//! zflate never reads or writes outside the buffers it is given: each input and output sits
//! against an inaccessible page, just before one and just after one, so that a stray access
//! faults at once instead of going unnoticed. Results are checked against libdeflate's.
//!
//! `ZFLATE_GUARD_TRIALS` sets the damaged variants per stream (default 40).

#![cfg(unix)]

#[path = "../benches/support/corpus.rs"]
mod corpus;

use libdeflater::{CompressionLvl, Compressor};
use zflate::{DecompressionError, Decompressor};

/// Where a [`Guarded`] buffer sits against its inaccessible page.
#[derive(Clone, Copy, Debug)]
enum Place {
    /// Its last byte is just before the page.
    End,
    /// Its first byte is just after the page.
    Start,
}

const PLACES: [Place; 2] = [Place::End, Place::Start];

/// A buffer between two inaccessible pages, against one of them.
struct Guarded {
    map: *mut u8,
    map_len: usize,
    data: *mut u8,
    len: usize,
}

impl Guarded {
    fn new(len: usize, place: Place) -> Self {
        // SAFETY: sysconf has no preconditions.
        let page = unsafe { libc::sysconf(libc::_SC_PAGESIZE) } as usize;
        let pages = len.div_ceil(page).max(1);
        let map_len = (pages + 2) * page;
        // SAFETY: an anonymous private mapping, checked below.
        let map = unsafe {
            libc::mmap(
                std::ptr::null_mut(),
                map_len,
                libc::PROT_READ | libc::PROT_WRITE,
                libc::MAP_PRIVATE | libc::MAP_ANONYMOUS,
                -1,
                0,
            )
        };
        assert_ne!(map, libc::MAP_FAILED, "mmap of {map_len} bytes");
        let map = map.cast::<u8>();
        // SAFETY: the first and last pages of the mapping, and the data between them.
        unsafe {
            assert_eq!(libc::mprotect(map.cast(), page, libc::PROT_NONE), 0);
            let last = map.add((pages + 1) * page);
            assert_eq!(libc::mprotect(last.cast(), page, libc::PROT_NONE), 0);
            let data = match place {
                Place::End => last.sub(len),
                Place::Start => map.add(page),
            };
            Self {
                map,
                map_len,
                data,
                len,
            }
        }
    }

    fn with(bytes: &[u8], place: Place) -> Self {
        let mut g = Self::new(bytes.len(), place);
        g.as_mut().copy_from_slice(bytes);
        g
    }

    fn as_ref(&self) -> &[u8] {
        // SAFETY: `len` readable bytes at `data`, owned by `self`.
        unsafe { std::slice::from_raw_parts(self.data, self.len) }
    }

    fn as_mut(&mut self) -> &mut [u8] {
        // SAFETY: `len` writable bytes at `data`, owned by `self`.
        unsafe { std::slice::from_raw_parts_mut(self.data, self.len) }
    }
}

impl Drop for Guarded {
    fn drop(&mut self) {
        // SAFETY: the mapping made in `new`.
        unsafe { libc::munmap(self.map.cast(), self.map_len) };
    }
}

#[derive(Clone, Copy, Debug)]
enum Format {
    Deflate,
    Zlib,
    Gzip,
}

/// zflate's result for `data` into a guarded output of `out_len` bytes, both at `place`.
fn zflate_guarded(
    d: &mut Decompressor,
    format: Format,
    data: &[u8],
    out_len: usize,
    place: Place,
) -> (Result<usize, DecompressionError>, Vec<u8>) {
    let input = Guarded::with(data, place);
    let mut out = Guarded::new(out_len, place);
    let r = match format {
        Format::Deflate => d.deflate_decompress(input.as_ref(), out.as_mut()),
        Format::Zlib => d.zlib_decompress(input.as_ref(), out.as_mut()),
        Format::Gzip => d.gzip_decompress(input.as_ref(), out.as_mut()),
    };
    let written = r.unwrap_or(0);
    (r, out.as_ref()[..written].to_vec())
}

/// libdeflate's result for the same, in ordinary buffers.
fn libdeflate(
    d: &mut libdeflater::Decompressor,
    format: Format,
    data: &[u8],
    out_len: usize,
) -> (Result<usize, DecompressionError>, Vec<u8>) {
    let mut out = vec![0; out_len];
    let r = match format {
        Format::Deflate => d.deflate_decompress(data, &mut out),
        Format::Zlib => d.zlib_decompress(data, &mut out),
        Format::Gzip => d.gzip_decompress(data, &mut out),
    }
    .map_err(|e| match e {
        libdeflater::DecompressionError::InsufficientSpace => DecompressionError::InsufficientSpace,
        libdeflater::DecompressionError::BadData => DecompressionError::BadData,
    });
    let written = r.unwrap_or(0);
    out.truncate(written);
    (r, out)
}

/// Checks one input at every placement and output size against libdeflate.
fn check(
    d: &mut Decompressor,
    lib: &mut libdeflater::Decompressor,
    format: Format,
    data: &[u8],
    len: usize,
    what: &str,
) {
    for out_len in [len, len.saturating_sub(1), len / 2, len + 300] {
        let want = libdeflate(lib, format, data, out_len);
        for place in PLACES {
            let got = zflate_guarded(d, format, data, out_len, place);
            assert_eq!(got, want, "{what}: {format:?}, output {out_len}, {place:?}");
        }
    }
}

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

#[test]
fn real_blocks_stay_in_bounds() {
    let mut d = Decompressor::new();
    let mut lib = libdeflater::Decompressor::new();
    for file in corpus::files() {
        eprintln!("{}: {} blocks", file.name, file.blocks.len());
        for (i, block) in file.blocks.iter().enumerate() {
            for place in PLACES {
                let (r, out) =
                    zflate_guarded(&mut d, Format::Deflate, &block.deflate, block.len, place);
                assert_eq!(r, Ok(block.len), "{} block {i} {place:?}", file.name);
                assert_eq!(
                    zflate::crc32(0, &out),
                    block.crc32,
                    "{} block {i}",
                    file.name
                );
            }
            // Short and long outputs too, on a sample of the blocks.
            if i % 16 == 0 {
                check(
                    &mut d,
                    &mut lib,
                    Format::Deflate,
                    &block.deflate,
                    block.len,
                    &format!("{} block {i}", file.name),
                );
            }
        }
    }
}

#[test]
fn damaged_streams_stay_in_bounds() {
    let trials: usize =
        std::env::var("ZFLATE_GUARD_TRIALS").map_or(40, |s| s.parse().expect("trials"));
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let mut d = Decompressor::new();
    let mut lib = libdeflater::Decompressor::new();
    let mut streams = 0;
    for len in [0, 1, 3, 100, 4096, 65_536] {
        let inputs: [Vec<u8>; 3] = [
            (0..len).map(|_| b"ACGT"[rng.below(4)]).collect(),
            (0..len).map(|_| rng.next() as u8).collect(),
            (0..len).map(|i| (i % 7) as u8).collect(),
        ];
        for data in &inputs {
            for level in [0, 1, 6, 12] {
                let mut c = Compressor::new(CompressionLvl::new(level).expect("level"));
                for format in [Format::Deflate, Format::Zlib, Format::Gzip] {
                    let mut z = vec![0; c.gzip_compress_bound(len) + 64];
                    let n = match format {
                        Format::Deflate => c.deflate_compress(data, &mut z),
                        Format::Zlib => c.zlib_compress(data, &mut z),
                        Format::Gzip => c.gzip_compress(data, &mut z),
                    }
                    .expect("compress");
                    z.truncate(n);
                    check(
                        &mut d,
                        &mut lib,
                        format,
                        &z,
                        len,
                        &format!("len {len} level {level}"),
                    );
                    streams += 1;
                    for trial in 0..trials {
                        let mut bad = z.clone();
                        match rng.below(3) {
                            0 if !bad.is_empty() => {
                                let i = rng.below(bad.len());
                                bad[i] ^= 1 << rng.below(8);
                            }
                            1 => bad.truncate(rng.below(bad.len() + 1)),
                            _ if !bad.is_empty() => {
                                let i = rng.below(bad.len());
                                bad[i] = rng.next() as u8;
                            }
                            _ => {}
                        }
                        check(
                            &mut d,
                            &mut lib,
                            format,
                            &bad,
                            len,
                            &format!("len {len} level {level} trial {trial}"),
                        );
                        streams += 1;
                    }
                }
            }
        }
    }
    eprintln!("{streams} streams, each at 4 output sizes and 2 placements");
}

/// Inputs for the compressor that put matches against the end: runs, short and long
/// periods, text, noise, and real blocks when the corpus is there.
fn compression_inputs() -> Vec<Vec<u8>> {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let mut inputs: Vec<Vec<u8>> = vec![vec![0; 70_000], vec![b'x'; 300]];
    for period in [1, 3, 4, 5, 9, 100, 258, 1000, 32_767] {
        let unit: Vec<u8> = (0..period).map(|_| rng.next() as u8).collect();
        for len in [period + 7, 40_000, 65_280] {
            inputs.push(unit.iter().copied().cycle().take(len).collect());
        }
    }
    inputs.push((0..65_280).map(|_| rng.next() as u8).collect());
    // Noise that turns into a repeat of itself near the end.
    let mut tail: Vec<u8> = (0..30_000).map(|_| rng.next() as u8).collect();
    tail.extend_from_within(29_000..30_000);
    inputs.push(tail);
    for file in corpus::files() {
        let mut d = Decompressor::new();
        for b in file.blocks.iter().step_by(97).take(8) {
            let mut out = vec![0; b.len];
            d.deflate_decompress(&b.deflate, &mut out)
                .expect("valid block");
            inputs.push(out);
        }
    }
    inputs
}

#[test]
fn compression_stays_in_bounds() {
    let inputs = compression_inputs();
    let mut d = Decompressor::new();
    for level in 0..=12 {
        let mut c = zflate::Compressor::new(zflate::CompressionLevel::new(level).expect("level"));
        for (i, data) in inputs.iter().enumerate() {
            let bound = zflate::deflate_compress_bound(data.len());
            // The size it compresses to, to try an output that fits exactly.
            let mut scratch = vec![0; bound];
            let n = c
                .deflate_compress(data, &mut scratch)
                .expect("fits in the bound");
            for place in PLACES {
                let input = Guarded::with(data, place);
                for (out_len, fits) in [(bound, true), (n, true), (n - 1, false)] {
                    let mut out = Guarded::new(out_len, Place::End);
                    let result = c.deflate_compress(input.as_ref(), out.as_mut());
                    let what = format!(
                        "level {level}, input {i} ({} bytes) at {place:?}, output {out_len}",
                        data.len()
                    );
                    if fits {
                        assert_eq!(result, Ok(n), "{what}");
                        let mut back = vec![0; data.len()];
                        assert_eq!(
                            d.deflate_decompress(&out.as_ref()[..n], &mut back),
                            Ok(data.len()),
                            "{what}"
                        );
                        assert!(back == *data, "{what}");
                    } else {
                        assert_eq!(
                            result,
                            Err(zflate::CompressionError::InsufficientSpace),
                            "{what}"
                        );
                    }
                }
            }
        }
    }
    eprintln!(
        "{} inputs, 13 levels, 2 input placements, 3 output sizes",
        inputs.len()
    );
}

/// `deflate_decompress_batch` with every input and output against guard pages, compared
/// with the single decoder in ordinary buffers: real blocks, and damaged streams, in
/// batches of 2 and 5.
#[test]
fn batches_stay_in_bounds() {
    let mut rng = Rng(0x5851_f42d_4c95_7f2d);
    let mut streams: Vec<(Vec<u8>, usize)> = Vec::new();
    for file in corpus::files() {
        for b in &file.blocks {
            streams.push((b.deflate.clone(), b.len));
            // One byte short, for some; not for the empty EOF block a whole file ends with.
            if streams.len().is_multiple_of(8) && b.len > 0 {
                streams.push((b.deflate.clone(), b.len - 1));
            }
        }
    }
    for len in [0, 100, 4096, 65_536] {
        let data: Vec<u8> = (0..len).map(|_| b"ACGT"[rng.below(4)]).collect();
        for level in [0, 1, 6, 12] {
            let mut c = Compressor::new(CompressionLvl::new(level).expect("level"));
            let mut z = vec![0; c.deflate_compress_bound(len)];
            let n = c.deflate_compress(&data, &mut z).expect("compress");
            z.truncate(n);
            for _ in 0..40 {
                let mut bad = z.clone();
                match rng.below(3) {
                    0 if !bad.is_empty() => {
                        let i = rng.below(bad.len());
                        bad[i] ^= 1 << rng.below(8);
                    }
                    1 => bad.truncate(rng.below(bad.len() + 1)),
                    _ if !bad.is_empty() => {
                        let i = rng.below(bad.len());
                        bad[i] = rng.next() as u8;
                    }
                    _ => {}
                }
                streams.push((bad, len));
            }
            streams.push((z, len));
        }
    }
    let mut single = Decompressor::new();
    let want: Vec<_> = streams
        .iter()
        .map(|(z, len)| {
            let mut out = vec![0; *len];
            (single.deflate_decompress(z, &mut out), out)
        })
        .collect();
    let mut d = Decompressor::new();
    for size in [2, 5] {
        for place in PLACES {
            for (k, group) in streams.chunks(size).enumerate() {
                let ins: Vec<Guarded> =
                    group.iter().map(|(z, _)| Guarded::with(z, place)).collect();
                let mut outs: Vec<Guarded> = group
                    .iter()
                    .map(|(_, len)| Guarded::new(*len, place))
                    .collect();
                let data: Vec<&[u8]> = ins.iter().map(Guarded::as_ref).collect();
                let mut out: Vec<&mut [u8]> = outs.iter_mut().map(Guarded::as_mut).collect();
                let mut got = vec![Ok(0); data.len()];
                d.deflate_decompress_batch(&data, &mut out, &mut got);
                for (j, r) in got.into_iter().enumerate() {
                    let i = k * size + j;
                    let (wr, wout) = &want[i];
                    assert_eq!(r, *wr, "stream {i}, batches of {size}, {place:?}");
                    if let Ok(n) = r {
                        assert!(
                            out[j][..n] == wout[..n],
                            "stream {i}, batches of {size}, {place:?}"
                        );
                    }
                }
            }
        }
    }
    eprintln!("{} streams in batches", streams.len());
}
