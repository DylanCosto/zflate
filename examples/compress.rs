//! Compresses a corpus file's BGZF blocks with zflate and libdeflate at each level, one core:
//! speed (the median of alternating rounds) and size, and checks zflate's output inflates
//! back. With `zflate` or `libdeflate` instead of `compare`, compresses with that one only,
//! for a profiler.
//!
//! `cargo run --release --example compress -- <compare|zflate|libdeflate> <file> [levels] [rounds]`
//!
//! `levels` is a comma-separated list (default `1,6,7`); `rounds` defaults to 5.

#[path = "../benches/support/corpus.rs"]
mod corpus;

use std::time::Instant;

use zflate::{CompressionLevel, Compressor, Decompressor, deflate_compress_bound};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let [_, which, path, rest @ ..] = &args[..] else {
        eprintln!("usage: compress <compare|zflate|libdeflate> <file> [levels] [rounds]");
        std::process::exit(2);
    };
    let levels: Vec<u8> = rest
        .first()
        .map_or("1,6,7", String::as_str)
        .split(',')
        .map(|l| l.parse().expect("level"))
        .collect();
    let rounds: usize = rest.get(1).map_or(5, |r| r.parse().expect("rounds"));
    // SAFETY: nothing else in this process reads or writes the environment.
    unsafe { std::env::set_var("ZFLATE_CORPUS", path) };
    let file = corpus::files().pop().expect("corpus file");
    let mut d = Decompressor::new();
    let blocks: Vec<Vec<u8>> = file
        .blocks
        .iter()
        .map(|b| {
            let mut out = vec![0; b.len];
            d.deflate_decompress(&b.deflate, &mut out)
                .expect("valid block");
            out
        })
        .collect();
    let bytes: usize = blocks.iter().map(Vec::len).sum();
    let mib = bytes as f64 / f64::from(1 << 20);
    println!("{} ({} blocks, {mib:.1} MiB)", file.name, blocks.len());
    let mut out = vec![0; deflate_compress_bound(65_536)];
    let mut check = vec![0; 65_536];

    for level in levels {
        let mut zflate = Compressor::new(CompressionLevel::new(level).expect("level 0 to 12"));
        let mut libdeflate = libdeflater::Compressor::new(
            libdeflater::CompressionLvl::new(i32::from(level)).expect("level 0 to 12"),
        );
        let mut zflate_once = |out: &mut [u8]| -> usize {
            blocks
                .iter()
                .map(|b| zflate.deflate_compress(b, out).expect("fits"))
                .sum()
        };
        let mut libdeflate_once = |out: &mut [u8]| -> usize {
            blocks
                .iter()
                .map(|b| libdeflate.deflate_compress(b, out).expect("fits"))
                .sum()
        };
        let time = |f: &mut dyn FnMut(&mut [u8]) -> usize, out: &mut [u8]| {
            let start = Instant::now();
            let size = f(out);
            (start.elapsed().as_secs_f64(), size)
        };
        match which.as_str() {
            "zflate" | "libdeflate" => {
                let f: &mut dyn FnMut(&mut [u8]) -> usize = if which == "zflate" {
                    &mut zflate_once
                } else {
                    &mut libdeflate_once
                };
                let (t, size) = (0..rounds).fold((0.0, 0), |(t, _), _| {
                    let (dt, size) = time(f, &mut out);
                    (t + dt, size)
                });
                println!(
                    "  level {level:2}: {which} {:8.1} MiB/s, {size} bytes",
                    mib * rounds as f64 / t
                );
            }
            "compare" => {
                let (mut tp, mut tl) = (Vec::new(), Vec::new());
                let (mut sp, mut sl) = (0, 0);
                for round in 0..rounds {
                    // Alternate which goes first.
                    for first in [round % 2 == 0, round % 2 != 0] {
                        if first {
                            let (t, s) = time(&mut zflate_once, &mut out);
                            tp.push(t);
                            sp = s;
                        } else {
                            let (t, s) = time(&mut libdeflate_once, &mut out);
                            tl.push(t);
                            sl = s;
                        }
                    }
                }
                // Check that every block inflates back.
                for b in &blocks {
                    let n = zflate.deflate_compress(b, &mut out).expect("fits");
                    let m = d
                        .deflate_decompress(&out[..n], &mut check)
                        .expect("inflates");
                    assert!(
                        check[..m] == b[..],
                        "level {level}: a block does not round-trip"
                    );
                }
                let median = |v: &mut Vec<f64>| {
                    v.sort_by(f64::total_cmp);
                    v[v.len() / 2]
                };
                let (rp, rl) = (mib / median(&mut tp), mib / median(&mut tl));
                println!(
                    "  level {level:2}: zflate {rp:8.1} MiB/s, libdeflate {rl:8.1} MiB/s, speed {:6.1}%;  size {sp} vs {sl}, {:7.3}%",
                    100.0 * rp / rl,
                    100.0 * sp as f64 / sl as f64
                );
            }
            _ => panic!("compare, zflate or libdeflate, not {which}"),
        }
    }
}
