//! zflate against libdeflate: the same answers for the same inputs, valid or not.

#[path = "../benches/support/corpus.rs"]
mod corpus;

use libdeflater::{CompressionLvl, Compressor};
use zflate::{DecompressionError, Decompressor};

/// A small deterministic xorshift generator, so failures reproduce.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// Inputs that exercise stored, fixed, and dynamic blocks, short and long matches, and
/// overlapping copies.
fn samples() -> Vec<(String, Vec<u8>)> {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let mut out = Vec::new();
    for len in [
        0, 1, 2, 3, 4, 7, 8, 9, 31, 257, 258, 259, 1000, 4096, 65_535, 65_536, 300_000,
    ] {
        out.push((format!("zeros/{len}"), vec![0; len]));
        out.push((
            format!("random/{len}"),
            (0..len).map(|_| rng.next() as u8).collect(),
        ));
        out.push((
            format!("acgt/{len}"),
            (0..len).map(|_| b"ACGT"[rng.below(4)]).collect(),
        ));
        let mut text = Vec::with_capacity(len);
        while text.len() < len {
            let word = [
                &b"chr1\t"[..],
                b"10583\t",
                b"rs58108140\t",
                b"G\tA\t",
                b"PASS\n",
            ];
            text.extend_from_slice(word[rng.below(word.len())]);
            if rng.below(16) == 0 {
                text.push(rng.next() as u8);
            }
        }
        text.truncate(len);
        out.push((format!("text/{len}"), text));
        // A run with a short period, for copies whose source overlaps their destination.
        out.push((
            format!("period/{len}"),
            (0..len).map(|i| (i % (1 + len % 7)) as u8).collect(),
        ));
    }
    out
}

#[derive(Clone, Copy, Debug)]
enum Format {
    Deflate,
    Zlib,
    Gzip,
}

const FORMATS: [Format; 3] = [Format::Deflate, Format::Zlib, Format::Gzip];

fn compress(c: &mut Compressor, format: Format, data: &[u8]) -> Vec<u8> {
    let bound = match format {
        Format::Deflate => c.deflate_compress_bound(data.len()),
        Format::Zlib => c.zlib_compress_bound(data.len()),
        Format::Gzip => c.gzip_compress_bound(data.len()),
    };
    let mut out = vec![0; bound];
    let n = match format {
        Format::Deflate => c.deflate_compress(data, &mut out),
        Format::Zlib => c.zlib_compress(data, &mut out),
        Format::Gzip => c.gzip_compress(data, &mut out),
    }
    .unwrap();
    out.truncate(n);
    out
}

fn zflate(
    d: &mut Decompressor,
    format: Format,
    data: &[u8],
    out: &mut [u8],
) -> Result<usize, DecompressionError> {
    match format {
        Format::Deflate => d.deflate_decompress(data, out),
        Format::Zlib => d.zlib_decompress(data, out),
        Format::Gzip => d.gzip_decompress(data, out),
    }
}

fn libdeflate(
    d: &mut libdeflater::Decompressor,
    format: Format,
    data: &[u8],
    out: &mut [u8],
) -> Result<usize, DecompressionError> {
    let r = match format {
        Format::Deflate => d.deflate_decompress(data, out),
        Format::Zlib => d.zlib_decompress(data, out),
        Format::Gzip => d.gzip_decompress(data, out),
    };
    r.map_err(|e| match e {
        libdeflater::DecompressionError::InsufficientSpace => DecompressionError::InsufficientSpace,
        libdeflater::DecompressionError::BadData => DecompressionError::BadData,
    })
}

#[test]
fn round_trips_every_level_and_format() {
    let mut d = Decompressor::new();
    for level in 0..=12 {
        let mut c = Compressor::new(CompressionLvl::new(level).unwrap());
        for (name, data) in samples() {
            for format in FORMATS {
                let z = compress(&mut c, format, &data);
                let ctx = format!("{name} level {level} {format:?}");

                let mut out = vec![0xaa; data.len()];
                assert_eq!(
                    zflate(&mut d, format, &z, &mut out),
                    Ok(data.len()),
                    "{ctx}"
                );
                assert!(out == data, "{ctx}: wrong output");

                let mut roomy = vec![0xaa; data.len() + 100];
                assert_eq!(
                    zflate(&mut d, format, &z, &mut roomy),
                    Ok(data.len()),
                    "{ctx}: roomy"
                );
                assert!(roomy[..data.len()] == data, "{ctx}: roomy: wrong output");

                if let Some(short) = data.len().checked_sub(1) {
                    let mut out = vec![0; short];
                    assert_eq!(
                        zflate(&mut d, format, &z, &mut out),
                        Err(DecompressionError::InsufficientSpace),
                        "{ctx}: one byte short"
                    );
                }
            }
        }
    }
}

/// Damaged variants per sample stream; `ZFLATE_DAMAGE_TRIALS` raises it for a longer run.
fn damage_trials() -> usize {
    std::env::var("ZFLATE_DAMAGE_TRIALS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(200)
}

/// Damaged streams must not panic, and must be accepted exactly when libdeflate accepts them,
/// with the same output.
#[test]
fn agrees_on_damaged_input() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let mut d = Decompressor::new();
    let mut l = libdeflater::Decompressor::new();
    let mut mismatches = Vec::new();
    for level in [0, 1, 6, 12] {
        let mut c = Compressor::new(CompressionLvl::new(level).unwrap());
        for (name, data) in samples().into_iter().filter(|(_, d)| d.len() <= 65_536) {
            for format in FORMATS {
                let z = compress(&mut c, format, &data);
                for trial in 0..damage_trials() {
                    let mut bad = z.clone();
                    match trial % 3 {
                        0 if !bad.is_empty() => {
                            let i = rng.below(bad.len());
                            bad[i] ^= 1 << rng.below(8);
                        }
                        1 => bad.truncate(rng.below(bad.len() + 1)),
                        _ => {
                            for _ in 0..1 + rng.below(4) {
                                if !bad.is_empty() {
                                    let i = rng.below(bad.len());
                                    bad[i] = rng.next() as u8;
                                }
                            }
                        }
                    }
                    let mut a = vec![0; data.len() + 64];
                    let mut b = vec![0; data.len() + 64];
                    let ra = zflate(&mut d, format, &bad, &mut a);
                    let rb = libdeflate(&mut l, format, &bad, &mut b);
                    let same = match (ra, rb) {
                        (Ok(n), Ok(m)) => n == m && a[..n] == b[..m],
                        (Err(_), Err(_)) => true,
                        _ => false,
                    };
                    if !same {
                        mismatches.push(format!(
                            "{name} level {level} {format:?} trial {trial}: zflate {ra:?}, libdeflate {rb:?}\n  \
                             good {}\n  bad  {}",
                            hex(&z),
                            hex(&bad),
                        ));
                    }
                }
            }
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} mismatches, e.g.:\n{}",
        mismatches.len(),
        mismatches[..mismatches.len().min(10)].join("\n")
    );
}

fn hex(data: &[u8]) -> String {
    let shown: String = data.iter().take(48).map(|b| format!("{b:02x}")).collect();
    if data.len() > 48 {
        format!("{shown}... ({} bytes)", data.len())
    } else {
        shown
    }
}

#[test]
fn checksums_match_libdeflate() {
    let mut rng = Rng(7);
    let data: Vec<u8> = (0..200_000).map(|_| rng.next() as u8).collect();
    for len in [
        0, 1, 7, 8, 15, 16, 63, 64, 65, 255, 256, 4097, 65_536, 200_000,
    ] {
        let s = &data[..len];
        assert_eq!(zflate::crc32(0, s), libdeflater::crc32(s), "crc32/{len}");
        assert_eq!(
            zflate::adler32(1, s),
            libdeflater::adler32(s),
            "adler32/{len}"
        );
        let (a, b) = s.split_at(len / 3);
        assert_eq!(
            zflate::crc32(zflate::crc32(0, a), b),
            libdeflater::crc32(s),
            "crc32 chained/{len}"
        );
        assert_eq!(
            zflate::adler32(zflate::adler32(1, a), b),
            libdeflater::adler32(s),
            "adler32 chained/{len}"
        );
    }
}

/// Each CRC32 path the CPU has, not just the one `crc32` picks, at every short length and
/// offset, and around every fold boundary of the long ones.
#[test]
fn every_crc32_implementation_matches_libdeflate() {
    let mut rng = Rng(11);
    let data: Vec<u8> = (0..(1 << 20) + 4096).map(|_| rng.next() as u8).collect();
    let mut lens: Vec<usize> = (0..=2048).collect();
    for boundary in [4096, 8192, 65_536, 1 << 20] {
        lens.extend(boundary - 70..boundary + 70);
    }
    let impls = zflate::crc32_implementations();
    eprintln!(
        "crc32 implementations: {:?}",
        impls.iter().map(|(n, _)| n).collect::<Vec<_>>()
    );
    for (name, f) in impls {
        let crc = |seed: u32, s: &[u8]| !f(!seed, s);
        for &len in &lens {
            for offset in if len <= 2048 { 0..64 } else { 0..3 } {
                let s = &data[offset..offset + len];
                assert_eq!(
                    crc(0, s),
                    libdeflater::crc32(s),
                    "{name} len {len} offset {offset}"
                );
            }
        }
        let s = &data[..100_000];
        for split in [0, 1, 15, 16, 17, 63, 64, 65, 511, 512, 513, 50_000] {
            let (a, b) = s.split_at(split);
            assert_eq!(
                crc(crc(0, a), b),
                libdeflater::crc32(s),
                "{name} chained at {split}"
            );
        }
    }
}

#[test]
fn inflates_real_bgzf_blocks() {
    let files = corpus::files();
    if files.is_empty() {
        eprintln!("no corpus files; set ZFLATE_CORPUS");
    }
    let mut d = Decompressor::new();
    for file in files {
        eprintln!("{}: {} blocks", file.name, file.blocks.len());
        let mut out = vec![0; 65_536];
        for (i, block) in file.blocks.iter().enumerate() {
            let n = d
                .deflate_decompress(&block.deflate, &mut out[..block.len])
                .unwrap_or_else(|e| panic!("{} block {i}: {e}", file.name));
            assert_eq!(n, block.len, "{} block {i}", file.name);
            assert_eq!(
                zflate::crc32(0, &out[..n]),
                block.crc32,
                "{} block {i}",
                file.name
            );
        }
    }
}

/// Around the fastloop's limits: every output size from just short of the data to well past
/// the room one iteration needs, and every cut of the input's last bytes. Each result must
/// be libdeflate's.
#[test]
fn agrees_at_buffer_edges() {
    let mut d = Decompressor::new();
    let mut l = libdeflater::Decompressor::new();
    let samples = samples();
    let picks = [
        "text/65536",
        "acgt/4096",
        "period/1000",
        "random/4096",
        "zeros/65536",
    ];
    for level in [1, 6, 12] {
        let mut c = Compressor::new(CompressionLvl::new(level).unwrap());
        for (name, data) in samples.iter().filter(|(n, _)| picks.contains(&n.as_str())) {
            let z = compress(&mut c, Format::Deflate, data);
            for out_len in data.len().saturating_sub(5)..data.len() + 320 {
                let (mut a, mut b) = (vec![0; out_len], vec![0; out_len]);
                let ra = zflate(&mut d, Format::Deflate, &z, &mut a);
                let rb = libdeflate(&mut l, Format::Deflate, &z, &mut b);
                assert_eq!(ra, rb, "{name} level {level} out {out_len}");
                if let Ok(n) = ra {
                    assert!(a[..n] == b[..n], "{name} level {level} out {out_len}");
                }
            }
            for cut in 1..64.min(z.len()) {
                let input = &z[..z.len() - cut];
                let (mut a, mut b) = (vec![0; data.len()], vec![0; data.len()]);
                let ra = zflate(&mut d, Format::Deflate, input, &mut a);
                let rb = libdeflate(&mut l, Format::Deflate, input, &mut b);
                assert_eq!(ra.is_ok(), rb.is_ok(), "{name} level {level} cut {cut}");
            }
        }
    }
}

/// `deflate_decompress_batch` against a call of `deflate_decompress` per stream: the same
/// result, error kind and output, for `streams` in batches of several sizes, so that lanes
/// take up new streams at every kind of point, after errors too.
fn check_batches(d: &mut Decompressor, streams: &[(Vec<u8>, usize)], what: &str) {
    let mut single = Decompressor::new();
    let want: Vec<_> = streams
        .iter()
        .map(|(z, len)| {
            let mut out = vec![0; *len];
            let r = single.deflate_decompress(z, &mut out);
            (r, out)
        })
        .collect();
    for size in [1, 2, 3, 7, streams.len().max(1)] {
        for (k, group) in streams.chunks(size).enumerate() {
            let data: Vec<&[u8]> = group.iter().map(|(z, _)| &z[..]).collect();
            let mut bufs: Vec<Vec<u8>> = group.iter().map(|(_, len)| vec![0; *len]).collect();
            let mut outs: Vec<&mut [u8]> = bufs.iter_mut().map(|b| &mut b[..]).collect();
            let mut got = vec![Ok(0); data.len()];
            d.deflate_decompress_batch(&data, &mut outs, &mut got);
            for (j, (r, out)) in got.into_iter().zip(&outs).enumerate() {
                let i = k * size + j;
                let (wr, wout) = &want[i];
                assert_eq!(r, *wr, "{what}: stream {i}, batches of {size}");
                if let Ok(n) = r {
                    assert!(
                        out[..n] == wout[..n],
                        "{what}: stream {i}, batches of {size}: output"
                    );
                }
            }
        }
    }
}

#[test]
fn batches_agree_with_single_streams() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let mut d = Decompressor::new();
    let mut streams = Vec::new();
    for level in [0, 1, 6, 12] {
        let mut c = Compressor::new(CompressionLvl::new(level).unwrap());
        for (_, data) in samples().into_iter().filter(|(_, d)| d.len() <= 65_536) {
            let z = compress(&mut c, Format::Deflate, &data);
            streams.push((z.clone(), data.len()));
            for trial in 0..damage_trials() / 4 {
                let mut bad = z.clone();
                match trial % 3 {
                    0 if !bad.is_empty() => {
                        let i = rng.below(bad.len());
                        bad[i] ^= 1 << rng.below(8);
                    }
                    1 => bad.truncate(rng.below(bad.len() + 1)),
                    _ => {
                        if !bad.is_empty() {
                            let i = rng.below(bad.len());
                            bad[i] = rng.next() as u8;
                        }
                    }
                }
                streams.push((bad, data.len() + rng.below(64)));
            }
        }
    }
    check_batches(&mut d, &streams, "samples");
}

#[test]
fn batches_agree_at_buffer_edges() {
    let mut d = Decompressor::new();
    let samples = samples();
    let picks = ["text/65536", "acgt/4096", "period/1000", "random/4096"];
    for level in [1, 6, 12] {
        let mut c = Compressor::new(CompressionLvl::new(level).unwrap());
        for (_, data) in samples.iter().filter(|(n, _)| picks.contains(&n.as_str())) {
            let z = compress(&mut c, Format::Deflate, data);
            let mut streams = Vec::new();
            for out_len in data.len().saturating_sub(5)..data.len() + 320 {
                streams.push((z.clone(), out_len));
            }
            for cut in 1..64.min(z.len()) {
                streams.push((z[..z.len() - cut].to_vec(), data.len()));
            }
            check_batches(&mut d, &streams, "edges");
        }
    }
}

#[test]
fn batches_inflate_real_bgzf_blocks() {
    let files = corpus::files();
    if files.is_empty() {
        eprintln!("no corpus files; set ZFLATE_CORPUS");
    }
    let mut d = Decompressor::new();
    for file in files {
        let streams: Vec<_> = file
            .blocks
            .iter()
            .map(|b| (b.deflate.clone(), b.len))
            .collect();
        check_batches(&mut d, &streams, &file.name);
    }
}
