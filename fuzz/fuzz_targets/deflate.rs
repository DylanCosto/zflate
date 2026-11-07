//! The compressor: any input, at any level and in any format, must compress within the
//! bound and inflate back to itself with zflate and with libdeflate; given exactly the room
//! it needs, it must write the same; one byte less must fail. Compressors are reused across
//! inputs, as callers reuse them; small inputs are also checked against a fresh one.
//!
//! The first byte picks the level (mod 13) and the format; the second, for some inputs,
//! grows the rest past 64 KiB, so that the window moves and the hash tables rebase.

#![no_main]

use std::cell::RefCell;

use libfuzzer_sys::fuzz_target;
use zflate::{
    CompressionLevel, Compressor, Decompressor, deflate_compress_bound, gzip_compress_bound,
    zlib_compress_bound,
};

thread_local! {
    static REUSED: RefCell<Vec<Compressor>> = RefCell::new(
        (0..=12)
            .map(|l| Compressor::new(CompressionLevel::new(l).expect("level")))
            .collect(),
    );
}

#[derive(Clone, Copy, Debug)]
enum Format {
    Deflate,
    Zlib,
    Gzip,
}

fn compress(c: &mut Compressor, format: Format, data: &[u8], out: &mut [u8]) -> Option<usize> {
    match format {
        Format::Deflate => c.deflate_compress(data, out),
        Format::Zlib => c.zlib_compress(data, out),
        Format::Gzip => c.gzip_compress(data, out),
    }
    .ok()
}

fn bound(format: Format, len: usize) -> usize {
    match format {
        Format::Deflate => deflate_compress_bound(len),
        Format::Zlib => zlib_compress_bound(len),
        Format::Gzip => gzip_compress_bound(len),
    }
}

/// `seed` repeated to `len` bytes, each copy changed a little, so that there are matches at
/// many offsets, near and far.
fn grow(seed: &[u8], len: usize) -> Vec<u8> {
    let mut data = Vec::with_capacity(len);
    let mut copy = 0u32;
    while data.len() < len {
        for (i, &b) in seed.iter().enumerate() {
            if data.len() == len {
                break;
            }
            let change = if (i as u32 + copy) % 61 == 0 {
                copy as u8
            } else {
                0
            };
            data.push(b ^ change);
        }
        copy += 1;
    }
    data
}

fuzz_target!(|input: &[u8]| {
    let [select, size, rest @ ..] = input else {
        return;
    };
    let level = usize::from(select % 13);
    let format = match (select / 13) % 3 {
        1 => Format::Zlib,
        2 => Format::Gzip,
        _ => Format::Deflate,
    };
    let grown;
    let data = if size & 0x80 != 0 && !rest.is_empty() {
        let lens = [
            65_535, 65_536, 65_537, 70_000, 98_303, 98_304, 131_073, 200_000,
        ];
        grown = grow(rest, lens[usize::from(size & 7)]);
        &grown[..]
    } else {
        rest
    };
    REUSED.with(|cs| {
        let c = &mut cs.borrow_mut()[level];
        let mut out = vec![0; bound(format, data.len())];
        let n = compress(c, format, data, &mut out).expect("fits in the bound");
        let out = &out[..n];

        // Inflate with both decoders.
        let mut back = vec![0; data.len()];
        let mut d = Decompressor::new();
        let ours = match format {
            Format::Deflate => d.deflate_decompress(out, &mut back),
            Format::Zlib => d.zlib_decompress(out, &mut back),
            Format::Gzip => d.gzip_decompress(out, &mut back),
        };
        assert_eq!(
            ours,
            Ok(data.len()),
            "zflate inflates level {level} {format:?}"
        );
        assert!(
            back == data,
            "zflate's inflation differs, level {level} {format:?}"
        );
        let mut ld = libdeflater::Decompressor::new();
        back.fill(0);
        let theirs = match format {
            Format::Deflate => ld.deflate_decompress(out, &mut back),
            Format::Zlib => ld.zlib_decompress(out, &mut back),
            Format::Gzip => ld.gzip_decompress(out, &mut back),
        };
        assert_eq!(
            theirs.ok(),
            Some(data.len()),
            "libdeflate inflates level {level} {format:?}"
        );
        assert!(
            back == data,
            "libdeflate's inflation differs, level {level} {format:?}"
        );

        // Exactly the room: the same output. One byte less: no output.
        let mut exact = vec![0; n];
        assert_eq!(
            compress(c, format, data, &mut exact),
            Some(n),
            "exact room, level {level}"
        );
        assert!(exact == out, "output depends on the room, level {level}");
        assert_eq!(
            compress(c, format, data, &mut exact[..n - 1]),
            None,
            "one byte short, level {level}"
        );

        // State left by earlier inputs must not change the output.
        if data.len() <= 4096 {
            let mut fresh = Compressor::new(CompressionLevel::new(level as u8).expect("level"));
            let mut again = vec![0; n];
            assert_eq!(compress(&mut fresh, format, data, &mut again), Some(n));
            assert!(
                again == out,
                "a reused compressor's output differs, level {level}"
            );
        }
    });
});
