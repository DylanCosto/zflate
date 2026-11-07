//! Inflates every BGZF block of a corpus file repeatedly, for a profiler.
//!
//! `cargo run --release --example profile -- <zflate|libdeflate|zflate+crc|libdeflate+crc> <file> [rounds]`

#[path = "../benches/support/corpus.rs"]
mod corpus;

use std::hint::black_box;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let [_, which, path, rest @ ..] = &args[..] else {
        eprintln!("usage: profile <zflate|libdeflate|zflate+crc|libdeflate+crc> <file> [rounds]");
        std::process::exit(2);
    };
    let rounds: usize = rest.first().map_or(20, |r| r.parse().expect("rounds"));
    // SAFETY: nothing else in this process reads or writes the environment.
    unsafe { std::env::set_var("ZFLATE_CORPUS", path) };
    let file = corpus::files().pop().expect("corpus file");
    let mut out = vec![0; 65_536];
    let mut zflate = zflate::Decompressor::new();
    let mut libdeflate = libdeflater::Decompressor::new();

    let start = Instant::now();
    if which == "zflate-batch" || which == "zflate-batch+crc" {
        // Blocks in batches of `BATCH` (default 32), as a BGZF reader might hand them over.
        let batch: usize = std::env::var("BATCH").map_or(32, |s| s.parse().expect("BATCH"));
        let mut bufs = vec![vec![0u8; 65_536]; batch];
        let mut results = vec![Ok(0); batch];
        for _ in 0..rounds {
            for group in file.blocks.chunks(batch) {
                let data: Vec<&[u8]> = group.iter().map(|b| &b.deflate[..]).collect();
                let mut outs: Vec<&mut [u8]> = bufs
                    .iter_mut()
                    .zip(group)
                    .map(|(buf, b)| &mut buf[..b.len])
                    .collect();
                let results = &mut results[..group.len()];
                zflate.deflate_decompress_batch(&data, &mut outs, results);
                for ((n, b), out) in results.iter().zip(group).zip(&outs) {
                    let n = n.unwrap();
                    if which.ends_with("+crc") {
                        assert_eq!(zflate::crc32(0, &out[..n]), b.crc32);
                    }
                    black_box(n);
                }
            }
        }
        let secs = start.elapsed().as_secs_f64();
        let bytes = (file.inflated_len() * rounds) as f64;
        println!(
            "{which}: {} blocks x {rounds}, {:.3} GiB/s",
            file.blocks.len(),
            bytes / secs / f64::from(1 << 30)
        );
        return;
    }
    for _ in 0..rounds {
        for b in &file.blocks {
            let n = match which.as_str() {
                "zflate+crc" => {
                    let n = zflate
                        .deflate_decompress(&b.deflate, &mut out[..b.len])
                        .unwrap();
                    assert_eq!(zflate::crc32(0, &out[..n]), b.crc32);
                    n
                }
                "libdeflate+crc" => {
                    let n = libdeflate
                        .deflate_decompress(&b.deflate, &mut out[..b.len])
                        .unwrap();
                    assert_eq!(libdeflater::crc32(&out[..n]), b.crc32);
                    n
                }
                "zflate" => zflate
                    .deflate_decompress(&b.deflate, &mut out[..b.len])
                    .unwrap(),
                "libdeflate" => libdeflate
                    .deflate_decompress(&b.deflate, &mut out[..b.len])
                    .unwrap(),
                other => panic!("unknown decoder {other}"),
            };
            black_box(n);
        }
    }
    let secs = start.elapsed().as_secs_f64();
    let bytes = (file.inflated_len() * rounds) as f64;
    println!(
        "{which}: {} blocks x {rounds}, {:.3} GiB/s",
        file.blocks.len(),
        bytes / secs / f64::from(1 << 30)
    );
}
