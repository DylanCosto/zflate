//! zflate against libdeflate (the target), zune-inflate (zflate's starting point), and zlib-rs
//! (the fastest other pure-Rust codec) on pasta's BGZF blocks.
//!
//! Throughput is in uncompressed bytes for every group, so the numbers compare across groups.

#[path = "support/corpus.rs"]
mod corpus;

use std::hint::black_box;
use std::time::Duration;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use flate2::{Compress, Compression, Decompress, FlushCompress, FlushDecompress};
use libdeflater::{CompressionLvl, Compressor};

/// The libdeflate levels pasta uses for zlib-style levels 1-9.
const LEVELS: [i32; 9] = [1, 2, 3, 5, 6, 7, 8, 10, 12];

/// Blocks per file to compress: compression is slow at the high levels, and BGZF blocks are
/// independent, so a sample is representative.
const COMPRESS_BLOCKS: usize = 32;

/// Each file's blocks, inflated.
fn inflated(file: &corpus::File) -> Vec<Vec<u8>> {
    let mut d = libdeflater::Decompressor::new();
    file.blocks
        .iter()
        .map(|b| {
            let mut out = vec![0; b.len];
            d.deflate_decompress(&b.deflate, &mut out).unwrap();
            out
        })
        .collect()
}

fn inflate(c: &mut Criterion) {
    for file in corpus::files() {
        let mut g = c.benchmark_group(format!("inflate/{}", file.name));
        g.throughput(Throughput::Bytes(file.inflated_len() as u64));
        g.sample_size(20);
        let mut out = vec![0; 65_536];

        let mut d = libdeflater::Decompressor::new();
        g.bench_function("libdeflate", |bench| {
            bench.iter(|| {
                for b in &file.blocks {
                    black_box(d.deflate_decompress(&b.deflate, &mut out[..b.len]).unwrap());
                }
            })
        });

        let mut d = zflate::Decompressor::new();
        g.bench_function("zflate", |bench| {
            bench.iter(|| {
                for b in &file.blocks {
                    black_box(d.deflate_decompress(&b.deflate, &mut out[..b.len]).unwrap());
                }
            })
        });

        g.bench_function("zune-inflate", |bench| {
            bench.iter(|| {
                for b in &file.blocks {
                    let options = zune_inflate::DeflateOptions::default().set_size_hint(b.len);
                    let mut d = zune_inflate::DeflateDecoder::new_with_options(&b.deflate, options);
                    black_box(d.decode_deflate().unwrap());
                }
            })
        });

        let mut d = Decompress::new(false);
        g.bench_function("zlib-rs", |bench| {
            bench.iter(|| {
                for b in &file.blocks {
                    d.reset(false);
                    d.decompress(&b.deflate, &mut out[..b.len], FlushDecompress::Finish)
                        .unwrap();
                    black_box(&out);
                }
            })
        });
        g.finish();
    }
}

fn checksums(c: &mut Criterion) {
    for file in corpus::files() {
        let blocks = inflated(&file);
        let mut g = c.benchmark_group(format!("crc32/{}", file.name));
        g.throughput(Throughput::Bytes(file.inflated_len() as u64));
        g.bench_function("libdeflate", |bench| {
            bench.iter(|| {
                blocks
                    .iter()
                    .map(|b| libdeflater::crc32(b))
                    .fold(0, u32::wrapping_add)
            })
        });
        g.bench_function("zflate", |bench| {
            bench.iter(|| {
                blocks
                    .iter()
                    .map(|b| zflate::crc32(0, b))
                    .fold(0, u32::wrapping_add)
            })
        });
        g.bench_function("crc32fast", |bench| {
            bench.iter(|| {
                blocks
                    .iter()
                    .map(|b| crc32fast::hash(b))
                    .fold(0, u32::wrapping_add)
            })
        });
        g.finish();

        let mut g = c.benchmark_group(format!("adler32/{}", file.name));
        g.throughput(Throughput::Bytes(file.inflated_len() as u64));
        g.bench_function("libdeflate", |bench| {
            bench.iter(|| {
                blocks
                    .iter()
                    .map(|b| libdeflater::adler32(b))
                    .fold(0, u32::wrapping_add)
            })
        });
        g.bench_function("zflate", |bench| {
            bench.iter(|| {
                blocks
                    .iter()
                    .map(|b| zflate::adler32(1, b))
                    .fold(0, u32::wrapping_add)
            })
        });
        g.finish();
    }
}

fn deflate(c: &mut Criterion) {
    for file in corpus::files() {
        let blocks: Vec<_> = inflated(&file).into_iter().take(COMPRESS_BLOCKS).collect();
        let len: usize = blocks.iter().map(Vec::len).sum();
        let mut out = vec![0; 2 * 65_536];
        for level in LEVELS {
            let mut g = c.benchmark_group(format!("deflate/{}/{level}", file.name));
            g.throughput(Throughput::Bytes(len as u64));
            g.sample_size(10);
            g.measurement_time(Duration::from_secs(3));

            let mut z = Compressor::new(CompressionLvl::new(level).unwrap());
            let size: usize = blocks
                .iter()
                .map(|b| z.deflate_compress(b, &mut out).unwrap())
                .sum();
            eprintln!(
                "deflate/{}/{level} libdeflate ratio {:.4}",
                file.name,
                size as f64 / len as f64
            );
            g.bench_function("libdeflate", |bench| {
                bench.iter(|| {
                    for b in &blocks {
                        black_box(z.deflate_compress(b, &mut out).unwrap());
                    }
                })
            });

            // zlib-rs has levels 1-9 only; it is a reference point, not the target.
            if level <= 9 {
                let mut z = Compress::new(Compression::new(level as u32), false);
                let mut run = |out: &mut [u8]| -> usize {
                    blocks
                        .iter()
                        .map(|b| {
                            z.reset();
                            z.compress(b, out, FlushCompress::Finish).unwrap();
                            z.total_out() as usize
                        })
                        .sum()
                };
                let size = run(&mut out);
                eprintln!(
                    "deflate/{}/{level} zlib-rs ratio {:.4}",
                    file.name,
                    size as f64 / len as f64
                );
                g.bench_function("zlib-rs", |bench| bench.iter(|| black_box(run(&mut out))));
            }
            g.finish();
        }
    }
}

criterion_group!(benches, inflate, checksums, deflate);
criterion_main!(benches);
