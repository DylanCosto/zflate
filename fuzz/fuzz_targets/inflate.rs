//! The decoder against libdeflate: any input, in any of the three formats, into an output
//! buffer of one of several sizes, must give libdeflate's result (the same output, or the
//! same error). A decompressor reused across inputs must answer as a fresh one does.
//!
//! The first byte picks the format (bits 3-4) and the output size (bits 0-2); the rest is
//! the compressed data.

#![no_main]

use std::cell::RefCell;

use libfuzzer_sys::fuzz_target;
use zflate::{DecompressionError, Decompressor};

thread_local! {
    static REUSED: RefCell<Decompressor> = RefCell::new(Decompressor::new());
}

#[derive(Clone, Copy)]
enum Format {
    Deflate,
    Zlib,
    Gzip,
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

fn libdeflate(format: Format, data: &[u8], out: &mut [u8]) -> Result<usize, DecompressionError> {
    let mut d = libdeflater::Decompressor::new();
    let r = match format {
        Format::Deflate => d.deflate_decompress(data, out),
        Format::Zlib => d.zlib_decompress(data, out),
        Format::Gzip => d.gzip_decompress(data, out),
    };
    r.map_err(|e| match e {
        libdeflater::DecompressionError::BadData => DecompressionError::BadData,
        libdeflater::DecompressionError::InsufficientSpace => DecompressionError::InsufficientSpace,
    })
}

fuzz_target!(|input: &[u8]| {
    let Some((&select, data)) = input.split_first() else {
        return;
    };
    let room = match select & 7 {
        0 => 0,
        1 => 1,
        2 => 17,
        3 => 300,
        4 => 4096,
        5 => 65_280,
        6 => 1 << 20,
        _ => data.len() * 4,
    };
    let format = match (select >> 3) & 3 {
        1 => Format::Zlib,
        2 => Format::Gzip,
        _ => Format::Deflate,
    };
    let mut ours = vec![0; room];
    let mut theirs = vec![0; room];
    let a = zflate(&mut Decompressor::new(), format, data, &mut ours);
    let b = libdeflate(format, data, &mut theirs);
    assert_eq!(a, b, "zflate against libdeflate");
    if let Ok(n) = a {
        assert!(ours[..n] == theirs[..n], "output differs from libdeflate's");
    }
    let mut again = vec![0; room];
    let c = REUSED.with(|d| zflate(&mut d.borrow_mut(), format, data, &mut again));
    assert_eq!(a, c, "a reused decompressor against a fresh one");
    if let Ok(n) = c {
        assert!(
            again[..n] == ours[..n],
            "a reused decompressor's output differs"
        );
    }
});
