//! The compressor at every level and in every format: its output must decompress to the
//! input with zflate, libdeflate and zlib, fit in the bound, and fail cleanly in a buffer one
//! byte too small. Inputs are the real blocks in `data/vectors.bin` and made-up edge cases.

use std::io::Read;

use zflate::{
    CompressionError, CompressionLevel, Compressor, Decompressor, deflate_compress_bound,
    gzip_compress_bound, zlib_compress_bound,
};

static VECTORS: &[u8] = include_bytes!("data/vectors.bin");

/// The inflated data of the test vectors (real BGZF blocks among them).
fn vector_inputs() -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    let mut rest = VECTORS;
    while let Some((_, r)) = rest.split_first() {
        let (zlen, r) = r.split_first_chunk::<4>().expect("compressed length");
        let (dlen, r) = r.split_first_chunk::<4>().expect("data length");
        let (_, r) = r.split_at(u32::from_le_bytes(*zlen) as usize);
        let (data, r) = r.split_at(u32::from_le_bytes(*dlen) as usize);
        out.push(data.to_vec());
        rest = r;
    }
    out
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

    fn bytes(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| self.next() as u8).collect()
    }
}

/// Text-like data: words from a small vocabulary, so matches of all lengths and offsets.
fn words(rng: &mut Rng, n: usize) -> Vec<u8> {
    const WORDS: [&str; 12] = [
        "the ", "chr20\t", "ACGT", "quality ", "read ", "\n", "60M ", "=\t", "GATTACA", "0/1:",
        "PASS;", "rs1234 ",
    ];
    let mut out = Vec::with_capacity(n + 16);
    while out.len() < n {
        out.extend_from_slice(WORDS[(rng.next() % WORDS.len() as u64) as usize].as_bytes());
    }
    out.truncate(n);
    out
}

/// Made-up inputs for the edge cases: sizes around the stored cutoffs, runs, periods around
/// the window size, incompressible data, and inputs long enough to move the window and to
/// split blocks.
fn synthetic_inputs() -> Vec<Vec<u8>> {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let mut v: Vec<Vec<u8>> = vec![Vec::new(), vec![7], b"abc".to_vec()];
    for n in [4, 5, 6, 7, 8, 9, 51, 52, 53, 56, 100, 258, 259, 1000] {
        v.push(words(&mut rng, n));
        v.push(rng.bytes(n));
    }
    v.push(vec![0; 100_000]);
    v.push(vec![b'A'; 65_280]);
    v.push(rng.bytes(65_280));
    v.push(rng.bytes(200_000));
    v.push(words(&mut rng, 65_280));
    v.push(words(&mut rng, 1_000_000));
    for period in [1, 2, 3, 4, 7, 258, 259, 300, 32_767, 32_768, 32_769, 40_000] {
        let unit = rng.bytes(period);
        v.push(unit.iter().copied().cycle().take(100_000).collect());
    }
    // Text, then noise, then text again: the block splitter should notice.
    let mut mixed = words(&mut rng, 30_000);
    mixed.extend(rng.bytes(30_000));
    mixed.extend(words(&mut rng, 30_000));
    mixed.extend(vec![0; 30_000]);
    v.push(mixed);
    v
}

#[derive(Clone, Copy, Debug)]
enum Format {
    Deflate,
    Zlib,
    Gzip,
}

fn compress(
    c: &mut Compressor,
    format: Format,
    data: &[u8],
    out: &mut [u8],
) -> Result<usize, CompressionError> {
    match format {
        Format::Deflate => c.deflate_compress(data, out),
        Format::Zlib => c.zlib_compress(data, out),
        Format::Gzip => c.gzip_compress(data, out),
    }
}

fn bound(format: Format, len: usize) -> usize {
    match format {
        Format::Deflate => deflate_compress_bound(len),
        Format::Zlib => zlib_compress_bound(len),
        Format::Gzip => gzip_compress_bound(len),
    }
}

/// Decompresses `z` with zflate, libdeflate, and (for zlib and gzip) zlib, checking each
/// gives `data`.
fn check_decodes(format: Format, z: &[u8], data: &[u8], what: &str) {
    let mut out = vec![0; data.len()];
    let mut d = Decompressor::new();
    let n = match format {
        Format::Deflate => d.deflate_decompress(z, &mut out),
        Format::Zlib => d.zlib_decompress(z, &mut out),
        Format::Gzip => d.gzip_decompress(z, &mut out),
    };
    assert_eq!(n, Ok(data.len()), "zflate, {what}");
    assert!(out == data, "zflate, {what}");

    let mut ld = libdeflater::Decompressor::new();
    out.fill(0);
    let n = match format {
        Format::Deflate => ld.deflate_decompress(z, &mut out),
        Format::Zlib => ld.zlib_decompress(z, &mut out),
        Format::Gzip => ld.gzip_decompress(z, &mut out),
    };
    assert_eq!(n.ok(), Some(data.len()), "libdeflate, {what}");
    assert!(out == data, "libdeflate, {what}");

    let mut zlib_out = Vec::new();
    match format {
        Format::Deflate => flate2::read::DeflateDecoder::new(z).read_to_end(&mut zlib_out),
        Format::Zlib => flate2::read::ZlibDecoder::new(z).read_to_end(&mut zlib_out),
        Format::Gzip => flate2::read::GzDecoder::new(z).read_to_end(&mut zlib_out),
    }
    .unwrap_or_else(|e| panic!("zlib, {what}: {e}"));
    assert!(zlib_out == data, "zlib, {what}");
}

fn round_trip(inputs: &[Vec<u8>], formats: &[Format]) {
    for level in 0..=12 {
        let mut c = Compressor::new(CompressionLevel::new(level).expect("level"));
        for (i, data) in inputs.iter().enumerate() {
            for &format in formats {
                let what = format!(
                    "level {level}, {format:?}, input {i} ({} bytes)",
                    data.len()
                );
                let mut out = vec![0; bound(format, data.len())];
                let n = compress(&mut c, format, data, &mut out)
                    .unwrap_or_else(|e| panic!("{what}: {e}"));
                check_decodes(format, &out[..n], data, &what);
                // One byte less must fail, and fail cleanly.
                let mut short = vec![0; n - 1];
                assert_eq!(
                    compress(&mut c, format, data, &mut short),
                    Err(CompressionError::InsufficientSpace),
                    "{what}, into {} bytes",
                    n - 1
                );
                // The output does not depend on the room it has.
                let mut exact = vec![0; n];
                assert_eq!(compress(&mut c, format, data, &mut exact), Ok(n), "{what}");
                assert!(exact == out[..n], "{what}");
            }
        }
    }
}

#[test]
fn synthetic_round_trips() {
    round_trip(
        &synthetic_inputs(),
        &[Format::Deflate, Format::Zlib, Format::Gzip],
    );
}

#[test]
fn real_block_round_trips() {
    // The vectors repeat the same blocks in three formats; once each is enough.
    let inputs = vector_inputs();
    let mut seen = std::collections::HashSet::new();
    let inputs: Vec<Vec<u8>> = inputs
        .into_iter()
        .filter(|d| seen.insert(d.clone()))
        .collect();
    round_trip(&inputs, &[Format::Deflate]);
}

/// A buffer too small for the first of many stored blocks fails cleanly, however many
/// follow. Level 0 once kept each failed block's header bits, and panicked when more than
/// 21 had failed and a short last block fit.
#[test]
fn stored_blocks_fail_cleanly_in_a_small_buffer() {
    let data = vec![0; 21 * usize::from(u16::MAX) + 10];
    let mut c = Compressor::new(CompressionLevel::new(0).expect("level"));
    for format in [Format::Deflate, Format::Zlib, Format::Gzip] {
        for len in [0, 1, 23, 100, 65_535, data.len()] {
            let mut out = vec![0; len];
            assert_eq!(
                compress(&mut c, format, &data, &mut out),
                Err(CompressionError::InsufficientSpace),
                "{format:?} into {len} bytes"
            );
        }
    }
}

/// The output is about as small as libdeflate's at each level: every input within 1% (and
/// a few bytes, for tiny ones), and all of them together within 0.2%. Checking only the
/// total once let a table error that cost 3% on VCF text through.
///
/// Level 5 is the exception: it searches less deeply than libdeflate's and takes no 3-byte
/// matches, for 23-26% more speed on real BAMs and 22-42% on other files, at -1% to +1.5%
/// size on them. Small synthetic inputs, like the repetitive text here, can lose more; it
/// gets 3%.
#[test]
fn sizes_track_libdeflate() {
    let mut inputs = vector_inputs();
    let mut seen = std::collections::HashSet::new();
    inputs.retain(|d| seen.insert(d.clone()));
    let mut rng = Rng(42);
    inputs.push(words(&mut rng, 1_000_000));
    for level in 0..=12 {
        let mut c = Compressor::new(CompressionLevel::new(level).expect("level"));
        let mut lc = libdeflater::Compressor::new(
            libdeflater::CompressionLvl::new(i32::from(level)).expect("level"),
        );
        let (mut ours, mut theirs) = (0, 0);
        let (per_input, total) = if level == 5 {
            (1.03, 1.03)
        } else {
            (1.01, 1.002)
        };
        for (i, data) in inputs.iter().enumerate() {
            let mut out = vec![0; deflate_compress_bound(data.len())];
            let n = c.deflate_compress(data, &mut out).expect("fits");
            let m = lc.deflate_compress(data, &mut out).expect("fits");
            assert!(
                n as f64 <= m as f64 * per_input + 8.0,
                "level {level}, input {i} ({} bytes): {n} bytes against libdeflate's {m}",
                data.len()
            );
            ours += n;
            theirs += m;
        }
        let ratio = ours as f64 / theirs as f64;
        eprintln!(
            "level {level}: {ours} bytes, libdeflate {theirs}, {:.3}%",
            100.0 * ratio
        );
        assert!(
            ratio < total,
            "level {level}: {ours} bytes against libdeflate's {theirs}"
        );
    }
}
