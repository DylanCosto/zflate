//! Small fixed inputs, compressed by zlib (`data/make_vectors.py`), decoded with no other
//! library involved, so that the test also runs under Miri: each stream whole, into exact,
//! short and long outputs, and damaged variants of it.
//!
//! `ZFLATE_VECTOR_TRIALS` sets the damaged variants per stream (default 8).

use zflate::{DecompressionError, Decompressor};

static VECTORS: &[u8] = include_bytes!("data/vectors.bin");

/// A stream: its format (0 raw DEFLATE, 1 zlib, 2 gzip), its bytes, and what they inflate to.
fn vectors() -> Vec<(u8, &'static [u8], &'static [u8])> {
    let mut out = Vec::new();
    let mut rest = VECTORS;
    while let Some((&format, r)) = rest.split_first() {
        let (zlen, r) = r.split_first_chunk::<4>().expect("compressed length");
        let (dlen, r) = r.split_first_chunk::<4>().expect("data length");
        let (z, r) = r.split_at(u32::from_le_bytes(*zlen) as usize);
        let (data, r) = r.split_at(u32::from_le_bytes(*dlen) as usize);
        out.push((format, z, data));
        rest = r;
    }
    out
}

fn decompress(
    d: &mut Decompressor,
    format: u8,
    z: &[u8],
    out: &mut [u8],
) -> Result<usize, DecompressionError> {
    match format {
        0 => d.deflate_decompress(z, out),
        1 => d.zlib_decompress(z, out),
        _ => d.gzip_decompress(z, out),
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
fn vectors_round_trip() {
    let mut d = Decompressor::new();
    let vectors = vectors();
    assert!(vectors.len() > 400, "{} vectors", vectors.len());
    for (i, &(format, z, data)) in vectors.iter().enumerate() {
        let mut out = vec![0; data.len() + 300];
        let n = decompress(&mut d, format, z, &mut out[..data.len()]);
        assert_eq!(n, Ok(data.len()), "vector {i}");
        assert!(out[..data.len()] == *data, "vector {i}");
        let n = decompress(&mut d, format, z, &mut out);
        assert_eq!(n, Ok(data.len()), "vector {i}, longer output");
        assert!(out[..data.len()] == *data, "vector {i}, longer output");
        if !data.is_empty() {
            let short = decompress(&mut d, format, z, &mut out[..data.len() - 1]);
            assert_eq!(
                short,
                Err(DecompressionError::InsufficientSpace),
                "vector {i}"
            );
        }
    }
}

#[test]
fn damaged_vectors_decode_safely() {
    let trials: usize =
        std::env::var("ZFLATE_VECTOR_TRIALS").map_or(8, |s| s.parse().expect("trials"));
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let mut d = Decompressor::new();
    for &(format, z, data) in &vectors() {
        // The real blocks are large for Miri: damage them less.
        let trials = if data.len() > 10_000 {
            trials / 4
        } else {
            trials
        };
        for _ in 0..trials {
            let mut bad = z.to_vec();
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
            for out_len in [data.len(), data.len().saturating_sub(1), data.len() + 300] {
                let mut out = vec![0; out_len];
                // Any result will do; it must only come without undefined behaviour.
                let _ = decompress(&mut d, format, &bad, &mut out);
            }
        }
    }
}

/// The raw DEFLATE vectors through `deflate_decompress_batch`, three at a time, against
/// `deflate_decompress`: whole, one byte short, longer, and damaged, so that Miri checks
/// the batch decoder's unsafe code too.
#[test]
fn vectors_decode_in_batches() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let mut streams: Vec<(Vec<u8>, usize)> = Vec::new();
    for &(format, z, data) in &vectors() {
        if format != 0 {
            continue;
        }
        streams.push((z.to_vec(), data.len()));
        streams.push((z.to_vec(), data.len().saturating_sub(1)));
        streams.push((z.to_vec(), data.len() + 300));
        let mut bad = z.to_vec();
        if !bad.is_empty() {
            let i = rng.below(bad.len());
            bad[i] ^= 1 << rng.below(8);
        }
        streams.push((bad, data.len()));
    }
    let (mut single, mut d) = (Decompressor::new(), Decompressor::new());
    for (k, group) in streams.chunks(3).enumerate() {
        let want: Vec<_> = group
            .iter()
            .map(|(z, len)| {
                let mut out = vec![0; *len];
                (single.deflate_decompress(z, &mut out), out)
            })
            .collect();
        let data: Vec<&[u8]> = group.iter().map(|(z, _)| &z[..]).collect();
        let mut bufs: Vec<Vec<u8>> = group.iter().map(|(_, len)| vec![0; *len]).collect();
        let mut outs: Vec<&mut [u8]> = bufs.iter_mut().map(|b| &mut b[..]).collect();
        let mut got = vec![Ok(0); data.len()];
        d.deflate_decompress_batch(&data, &mut outs, &mut got);
        for (j, (r, (wr, wout))) in got.into_iter().zip(&want).enumerate() {
            assert_eq!(r, *wr, "stream {}", 3 * k + j);
            if let Ok(n) = r {
                assert!(outs[j][..n] == wout[..n], "stream {}", 3 * k + j);
            }
        }
    }
}

/// Compresses small inputs at every level and inflates them back, with zflate alone, so
/// that Miri checks the compressor's unsafe code too: some real blocks, a run, a pattern
/// that repeats to the very end, and an input over 64 KiB, which moves the window. Under
/// Miri, a representative set of levels, to finish in reasonable time.
#[test]
fn compressed_vectors_round_trip() {
    use zflate::{CompressionLevel, Compressor, deflate_compress_bound};

    let mut inputs: Vec<Vec<u8>> = vectors()
        .into_iter()
        .filter(|(_, _, data)| (200..4000).contains(&data.len()))
        .take(3)
        .map(|(_, _, data)| data.to_vec())
        .collect();
    inputs.push(vec![0; 40_000]);
    inputs.push(b"GATTACA".iter().copied().cycle().take(2000).collect());
    let levels: &[u8] = if cfg!(miri) {
        &[0, 1, 2, 5, 8, 10]
    } else {
        &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    };
    let mut d = Decompressor::new();
    let mut check = |level: u8, i: usize, data: &[u8]| {
        let mut c = Compressor::new(CompressionLevel::new(level).expect("level"));
        let mut out = vec![0; deflate_compress_bound(data.len())];
        let n = c
            .deflate_compress(data, &mut out)
            .expect("fits in the bound");
        let mut back = vec![0; data.len()];
        assert_eq!(
            d.deflate_decompress(&out[..n], &mut back),
            Ok(data.len()),
            "level {level}, input {i}"
        );
        assert!(back == *data, "level {level}, input {i}");
    };
    for &level in levels {
        for (i, data) in inputs.iter().enumerate() {
            check(level, i, data);
        }
    }
    // Past 64 KiB, the hash chains' window moves, and level 1 uses its rolling table,
    // which shares the short table's allocation.
    let large: Vec<u8> = b"ACGTTGCAAC".iter().copied().cycle().take(70_000).collect();
    for level in [1, 2, 5] {
        check(level, inputs.len(), &large);
    }
}
