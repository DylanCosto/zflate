//! Adler-32 alone, across buffer sizes: short buffers measure setup and tail handling, long ones
//! the main loop.

use std::hint::black_box;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

const SIZES: [usize; 13] = [
    4,
    16,
    32,
    64,
    128,
    256,
    512,
    1024,
    2048,
    4096,
    16_384,
    65_536,
    1 << 20,
];

fn adler32(c: &mut Criterion) {
    let data: Vec<u8> = (0..SIZES[SIZES.len() - 1])
        .map(|i| (i as u32).wrapping_mul(2_654_435_761).to_le_bytes()[3])
        .collect();
    let mut g = c.benchmark_group("adler32");
    g.warm_up_time(Duration::from_millis(500));
    g.measurement_time(Duration::from_secs(2));
    for len in SIZES {
        let s = &data[..len];
        g.throughput(Throughput::Bytes(len as u64));
        g.bench_with_input(BenchmarkId::new("libdeflate", len), s, |b, s| {
            b.iter(|| libdeflater::adler32(black_box(s)))
        });
        g.bench_with_input(BenchmarkId::new("zflate", len), s, |b, s| {
            b.iter(|| zflate::adler32(1, black_box(s)))
        });
        for (name, f) in zflate::adler32_implementations() {
            g.bench_with_input(
                BenchmarkId::new(format!("zflate-{name}"), len),
                s,
                |b, s| b.iter(|| f(1, black_box(s))),
            );
        }
    }
    g.finish();
}

criterion_group!(benches, adler32);
criterion_main!(benches);
