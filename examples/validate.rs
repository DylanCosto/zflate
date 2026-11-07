//! Checks zflate on whole real files: every BGZF block of a BAM, BCF or bgzipped file against
//! its CRC32 and length, or every gzip block of a CRAM 3 file against libdeflate.
//!
//! `cargo run --release --example validate -- <bgzf|cram>[-time] <file> [threads]`
//!
//! A block that fails is reported with its file offset, and with whether libdeflate decodes
//! it, which tells a bad file from a zflate bug.
//!
//! With `-time`, each worker also times both decoders on every chunk of blocks it checked,
//! one after the other on the same core, alternating which goes first: decompression and
//! the CRC32 check, as pasta reads BGZF (its own CRC32 function for each) and CRAM's gzip
//! blocks (whose CRC32 the gzip decoders check themselves).

use std::fs::File;
use std::io::{self, BufReader, Read};
use std::sync::mpsc;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

/// A run of whole blocks read from the file, and where each starts in it.
struct Job {
    /// Where `data` starts in the file.
    offset: u64,
    data: Vec<u8>,
    /// Each block's range in `data`.
    blocks: Vec<(usize, usize)>,
}

/// What the workers found.
#[derive(Default)]
struct Tally {
    blocks: u64,
    compressed: u64,
    inflated: u64,
    /// CRAM: blocks by compression method.
    methods: [u64; 16],
    failures: Vec<String>,
    /// With timing: the bytes timed, and each decoder's time on them.
    timed: u64,
    zflate_ns: u128,
    libdeflate_ns: u128,
}

impl Tally {
    fn add(&mut self, other: Tally) {
        self.blocks += other.blocks;
        self.compressed += other.compressed;
        self.inflated += other.inflated;
        for (a, b) in self.methods.iter_mut().zip(other.methods) {
            *a += b;
        }
        self.failures.extend(other.failures);
        self.timed += other.timed;
        self.zflate_ns += other.zflate_ns;
        self.libdeflate_ns += other.libdeflate_ns;
    }
}

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let [_, mode, path, rest @ ..] = &args[..] else {
        eprintln!("usage: validate <bgzf|cram> <file> [threads]");
        std::process::exit(2);
    };
    let threads: usize = rest.first().map_or(4, |t| t.parse().expect("threads"));
    let file = BufReader::with_capacity(1 << 20, File::open(path)?);
    let start = Instant::now();
    let (tx, rx) = mpsc::sync_channel::<Job>(2 * threads);
    let rx = Arc::new(Mutex::new(rx));
    let (cram, time) = match mode.as_str() {
        "bgzf" => (false, false),
        "cram" => (true, false),
        "bgzf-time" => (false, true),
        "cram-time" => (true, true),
        other => panic!("unknown mode {other}"),
    };
    let workers: Vec<_> = (0..threads)
        .map(|_| {
            let rx = Arc::clone(&rx);
            thread::spawn(move || {
                let mut tally = Tally::default();
                let mut zflate = zflate::Decompressor::new();
                let mut libdeflate = libdeflater::Decompressor::new();
                let mut out = vec![0; 1 << 21];
                let mut jobs = 0u64;
                loop {
                    let job = rx.lock().expect("lock poisoned").recv();
                    let Ok(job) = job else { break };
                    jobs += 1;
                    for &(at, len) in &job.blocks {
                        let block = &job.data[at..at + len];
                        let offset = job.offset + at as u64;
                        if cram {
                            check_cram_block(
                                block,
                                offset,
                                &mut zflate,
                                &mut libdeflate,
                                &mut tally,
                            );
                        } else {
                            check_bgzf_block(
                                block,
                                offset,
                                &mut zflate,
                                &mut libdeflate,
                                &mut tally,
                            );
                        }
                    }
                    if time {
                        let zflate_first = jobs.is_multiple_of(2);
                        for run in 0..2 {
                            let use_zflate = (run == 0) == zflate_first;
                            let start = Instant::now();
                            let bytes = time_job(
                                &job,
                                cram,
                                use_zflate,
                                &mut zflate,
                                &mut libdeflate,
                                &mut out,
                            );
                            let ns = start.elapsed().as_nanos();
                            if use_zflate {
                                tally.zflate_ns += ns;
                                tally.timed += bytes;
                            } else {
                                tally.libdeflate_ns += ns;
                            }
                        }
                    }
                }
                tally
            })
        })
        .collect();

    let read = if cram {
        read_cram(file, &tx)
    } else {
        read_bgzf(file, &tx)
    };
    drop(tx);
    let mut tally = Tally::default();
    for w in workers {
        tally.add(w.join().expect("worker panicked"));
    }
    read?;
    let secs = start.elapsed().as_secs_f64();
    println!(
        "{path}: {} blocks, {:.2} GB compressed, {:.2} GB inflated, {} failures, {:.0} s",
        tally.blocks,
        tally.compressed as f64 / 1e9,
        tally.inflated as f64 / 1e9,
        tally.failures.len(),
        secs
    );
    if cram {
        let names = [
            "raw", "gzip", "bzip2", "lzma", "rans4x8", "rans4x16", "arith", "fqzcomp", "tok3",
        ];
        let counts: Vec<String> = tally
            .methods
            .iter()
            .enumerate()
            .filter(|&(_, &n)| n > 0)
            .map(|(m, n)| format!("{} {n}", names.get(m).unwrap_or(&"other")))
            .collect();
        println!("  CRAM blocks by method: {}", counts.join(", "));
    }
    if time {
        let rate = |ns: u128| tally.timed as f64 / (ns as f64 / 1e9) / f64::from(1 << 30);
        let (p, l) = (rate(tally.zflate_ns), rate(tally.libdeflate_ns));
        println!(
            "  per core, with CRC32: zflate {p:.3} GiB/s, libdeflate {l:.3} GiB/s: {:.1}% ({:.2} GB timed)",
            100.0 * p / l,
            tally.timed as f64 / 1e9
        );
    }
    for f in tally.failures.iter().take(20) {
        println!("  FAIL {f}");
    }
    if !tally.failures.is_empty() {
        std::process::exit(1);
    }
    Ok(())
}

/// Decodes a job's blocks with one decoder, checking each CRC32 as pasta does, into `out`;
/// returns the bytes decoded. The blocks were already checked.
fn time_job(
    job: &Job,
    cram: bool,
    use_zflate: bool,
    zflate: &mut zflate::Decompressor,
    libdeflate: &mut libdeflater::Decompressor,
    out: &mut Vec<u8>,
) -> u64 {
    let mut bytes = 0;
    for &(at, len) in &job.blocks {
        let block = &job.data[at..at + len];
        if cram {
            let (method, size, raw, start) = cram_block_fields(block).expect("parsed when read");
            if method != 1 {
                continue;
            }
            if out.len() < raw {
                out.resize(raw, 0);
            }
            let gz = &block[start..start + size];
            let n = if use_zflate {
                zflate
                    .gzip_decompress(gz, &mut out[..raw])
                    .expect("checked")
            } else {
                libdeflate
                    .gzip_decompress(gz, &mut out[..raw])
                    .expect("checked")
            };
            bytes += n as u64;
        } else {
            let xlen = usize::from(u16::from_le_bytes([block[10], block[11]]));
            let trailer = &block[block.len() - 8..];
            let crc = u32::from_le_bytes(trailer[..4].try_into().expect("4 bytes"));
            let isize = u32::from_le_bytes(trailer[4..].try_into().expect("4 bytes")) as usize;
            let deflate = &block[12 + xlen..block.len() - 8];
            let dst = &mut out[..isize];
            let ok = if use_zflate {
                zflate.deflate_decompress(deflate, dst).expect("checked");
                zflate::crc32(0, dst) == crc
            } else {
                libdeflate
                    .deflate_decompress(deflate, dst)
                    .expect("checked");
                libdeflater::crc32(dst) == crc
            };
            assert!(ok, "checked");
            bytes += isize as u64;
        }
    }
    std::hint::black_box(&out);
    bytes
}

/// Splits a BGZF file into jobs of whole blocks, about 16 MB each.
fn read_bgzf(mut file: impl Read, tx: &mpsc::SyncSender<Job>) -> io::Result<()> {
    const CHUNK: usize = 16 << 20;
    let mut offset = 0u64;
    let mut carry: Vec<u8> = Vec::new();
    loop {
        let mut data = std::mem::take(&mut carry);
        let have = data.len();
        data.resize(have + CHUNK, 0);
        let n = read_full(&mut file, &mut data[have..])?;
        data.truncate(have + n);
        if data.is_empty() {
            return Ok(());
        }
        let mut blocks = Vec::new();
        let mut at = 0;
        while let Some(size) = bgzf_block_size(&data[at..]) {
            if at + size > data.len() {
                break;
            }
            blocks.push((at, size));
            at += size;
        }
        if n == 0 && at < data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("not a BGZF block at offset {}", offset + at as u64),
            ));
        }
        carry = data[at..].to_vec();
        data.truncate(at);
        let job = Job {
            offset,
            data,
            blocks,
        };
        offset += at as u64;
        if tx.send(job).is_err() {
            return Ok(());
        }
    }
}

/// The size of the BGZF block that starts `data`, if its header is all there.
fn bgzf_block_size(data: &[u8]) -> Option<usize> {
    if data.len() < 18 || data[..4] != [0x1f, 0x8b, 8, 4] {
        return None;
    }
    // The BC subfield is the first (and in practice only) extra subfield.
    if data[12..14] != *b"BC" {
        return None;
    }
    Some(usize::from(u16::from_le_bytes([data[16], data[17]])) + 1)
}

fn check_bgzf_block(
    block: &[u8],
    offset: u64,
    zflate: &mut zflate::Decompressor,
    libdeflate: &mut libdeflater::Decompressor,
    tally: &mut Tally,
) {
    let xlen = usize::from(u16::from_le_bytes([block[10], block[11]]));
    let trailer = &block[block.len() - 8..];
    let crc = u32::from_le_bytes(trailer[..4].try_into().expect("4 bytes"));
    let isize = u32::from_le_bytes(trailer[4..].try_into().expect("4 bytes")) as usize;
    let deflate = &block[12 + xlen..block.len() - 8];
    let mut out = vec![0; isize];
    let r = zflate.deflate_decompress(deflate, &mut out);
    tally.blocks += 1;
    tally.compressed += block.len() as u64;
    tally.inflated += isize as u64;
    if r != Ok(isize) || zflate::crc32(0, &out) != crc {
        let mut want = vec![0; isize];
        let lib = libdeflate.deflate_decompress(deflate, &mut want);
        let same = lib.is_ok() && want == out;
        tally.failures.push(format!(
            "BGZF block at offset {offset}: zflate {r:?}, libdeflate {lib:?}, outputs {}",
            if same { "equal" } else { "differ" }
        ));
    }
}

/// Reads a CRAM 3 file's containers and sends each one's blocks.
fn read_cram(mut file: impl Read, tx: &mpsc::SyncSender<Job>) -> io::Result<()> {
    let mut definition = [0u8; 26];
    file.read_exact(&mut definition)?;
    if definition[..4] != *b"CRAM" || definition[4] != 3 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("not CRAM 3: {:?}", &definition[..6]),
        ));
    }
    let mut offset = 26u64;
    loop {
        let mut length = [0u8; 4];
        match file.read_exact(&mut length) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(()),
            Err(e) => return Err(e),
        }
        let length = i32::from_le_bytes(length) as usize;
        let mut header = Reader::new(&mut file);
        header.itf8()?; // reference sequence id
        header.itf8()?; // alignment start
        header.itf8()?; // alignment span
        header.itf8()?; // records
        header.ltf8()?; // record counter
        header.ltf8()?; // bases
        let num_blocks = header.itf8()? as usize;
        let landmarks = header.itf8()?;
        for _ in 0..landmarks {
            header.itf8()?;
        }
        header.bytes(4)?; // CRC32
        let header_len = 4 + header.read as u64;
        let mut data = vec![0; length];
        file.read_exact(&mut data)?;
        let data_offset = offset + header_len;
        offset = data_offset + length as u64;
        let mut blocks = Vec::with_capacity(num_blocks);
        let mut at = 0;
        for _ in 0..num_blocks {
            let len = cram_block_len(&data[at..]).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("bad CRAM block at offset {}", data_offset + at as u64),
                )
            })?;
            blocks.push((at, len));
            at += len;
        }
        let job = Job {
            offset: data_offset,
            data,
            blocks,
        };
        if tx.send(job).is_err() {
            return Ok(());
        }
    }
}

/// A CRAM block's fields: method, compressed size and raw size, and where its data starts.
fn cram_block_fields(block: &[u8]) -> Option<(u8, usize, usize, usize)> {
    let method = *block.first()?;
    let mut r = Reader::new(block.get(2..)?);
    r.itf8().ok()?; // content id
    let size = r.itf8().ok()? as usize;
    let raw = r.itf8().ok()? as usize;
    Some((method, size, raw, 2 + r.read))
}

/// The length of the CRAM 3 block that starts `data`, CRC32 included.
fn cram_block_len(data: &[u8]) -> Option<usize> {
    let (_, size, _, start) = cram_block_fields(data)?;
    let len = start + size + 4;
    (len <= data.len()).then_some(len)
}

fn check_cram_block(
    block: &[u8],
    offset: u64,
    zflate: &mut zflate::Decompressor,
    libdeflate: &mut libdeflater::Decompressor,
    tally: &mut Tally,
) {
    let (method, size, raw, start) = cram_block_fields(block).expect("parsed when read");
    tally.blocks += 1;
    tally.methods[usize::from(method & 15)] += 1;
    let crc = u32::from_le_bytes(block[block.len() - 4..].try_into().expect("4 bytes"));
    if zflate::crc32(0, &block[..block.len() - 4]) != crc {
        tally.failures.push(format!(
            "CRAM block at offset {offset}: block CRC32 mismatch"
        ));
        return;
    }
    if method != 1 {
        return;
    }
    let gz = &block[start..start + size];
    tally.compressed += size as u64;
    tally.inflated += raw as u64;
    let mut out = vec![0; raw];
    let mut want = vec![0; raw];
    let r = zflate.gzip_decompress(gz, &mut out);
    let lib = libdeflate.gzip_decompress(gz, &mut want);
    let agree = match (&r, &lib) {
        (Ok(a), Ok(b)) => a == b && out == want,
        (Err(_), Err(_)) => true,
        _ => false,
    };
    if !agree || r != Ok(raw) {
        tally.failures.push(format!(
            "CRAM gzip block at offset {offset}: zflate {r:?}, libdeflate {lib:?}, raw size {raw}"
        ));
    }
}

/// Reads ITF8 and LTF8 integers (CRAM 3 section 2.3), counting the bytes read.
struct Reader<R> {
    inner: R,
    read: usize,
}

impl<R: Read> Reader<R> {
    fn new(inner: R) -> Self {
        Self { inner, read: 0 }
    }

    fn byte(&mut self) -> io::Result<u8> {
        let mut b = [0u8; 1];
        self.inner.read_exact(&mut b)?;
        self.read += 1;
        Ok(b[0])
    }

    fn bytes(&mut self, n: usize) -> io::Result<()> {
        for _ in 0..n {
            self.byte()?;
        }
        Ok(())
    }

    fn itf8(&mut self) -> io::Result<i32> {
        let b0 = u32::from(self.byte()?);
        let extra = match b0 {
            0..0x80 => return Ok(b0 as i32),
            0x80..0xc0 => 1,
            0xc0..0xe0 => 2,
            0xe0..0xf0 => 3,
            _ => 4,
        };
        // Four extra bytes keep four bits of the first, not three.
        let mut v = b0 & (0xff >> (extra + 1).min(4));
        for i in 0..extra {
            let b = u32::from(self.byte()?);
            v = if extra == 4 && i == 3 {
                (v << 4) | (b & 0x0f)
            } else {
                (v << 8) | b
            };
        }
        Ok(v as i32)
    }

    fn ltf8(&mut self) -> io::Result<i64> {
        let b0 = self.byte()?;
        let extra = b0.leading_ones() as usize;
        let mut v = if extra >= 8 {
            0
        } else {
            u64::from(b0 & (0xff >> (extra + 1)))
        };
        for _ in 0..extra {
            v = (v << 8) | u64::from(self.byte()?);
        }
        Ok(v as i64)
    }
}

/// Reads until `buf` is full or the input ends; returns the bytes read.
fn read_full(file: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match file.read(&mut buf[n..])? {
            0 => break,
            k => n += k,
        }
    }
    Ok(n)
}
