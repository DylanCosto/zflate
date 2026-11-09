# zflate

A pure-Rust DEFLATE (RFC 1951), zlib (RFC 1950) and gzip (RFC 1952) codec, ported from
[libdeflate](https://github.com/ebiggers/libdeflate) and tuned for genomics data: the BGZF
blocks of BAM, BCF and bgzipped VCF files, and CRAM's gzip blocks.

- **libdeflate's API**: whole buffers in, whole buffers out, no streaming. That is all BGZF
  and CRAM need, and it spares the decoder a streaming inflater's bookkeeping.
- **Decompression faster than libdeflate's** on every BGZF file tested, on x86-64 and Apple
  silicon.
- **Compression at libdeflate's sizes**, at every level from 0 to 12, with level 5 tuned
  for speed.
- **No C code and no dependencies.** CPU features are detected at run time.

## Use

```rust
use zflate::{CompressionLevel, Compressor, Decompressor, deflate_compress_bound};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let data = b"ACGTTGCA".repeat(1000);

    let mut compressor = Compressor::new(CompressionLevel::new(5).expect("0 to 12"));
    let mut compressed = vec![0; deflate_compress_bound(data.len())];
    let len = compressor.deflate_compress(&data, &mut compressed)?;
    compressed.truncate(len);

    let mut decompressor = Decompressor::new();
    let mut out = vec![0; data.len()];
    let len = decompressor.deflate_decompress(&compressed, &mut out)?;
    assert_eq!(out[..len], data);
    Ok(())
}
```

- `zlib_compress` / `zlib_decompress` and `gzip_compress` / `gzip_decompress` wrap the
  same codecs; `crc32` and `adler32` are exported too.
- A buffer of the `*_compress_bound` length always fits. A smaller one may fail with
  `CompressionError::InsufficientSpace`.
- The decompressor writes into a buffer you size: BGZF's and gzip's trailers give the
  length.
- Keep one `Compressor` per level and one `Decompressor` per thread, and reuse them: their
  tables are allocated once. Both are `Send` and `Sync`. `try_new` returns an error rather
  than aborting if the allocation fails.
- `Decompressor::deflate_decompress_batch` takes several independent streams, such as a
  file's BGZF blocks. On 64-bit Arm it decodes two at a time.

## Levels

| Level | Method | Size against libdeflate's |
|---|---|---|
| 0 | stored blocks | the same |
| 1 | hash table, greedy | within 0.1%, slightly smaller |
| 2–4 | hash chains, greedy | within 0.05% |
| 5 | hash chains, lazy, tuned for speed | −1% to +1.5% on BAM, BCF and VCF |
| 6–9 | hash chains, lazy | within 0.05% |
| 10–12 | near-optimal parsing | within 0.05% |

Level 5 searches less deeply than libdeflate's and takes no 3-byte matches: on BAM, BCF and
VCF it is 1.2–1.8 times as fast as libdeflate's level 5, at about the same size. Other data
can grow more: CRAM's gzip blocks from unaligned HiFi reads came out 8.8% larger than
libdeflate's level 5 ([benchmarks](docs/LONG-READ-CODECS.md)).

## Speed

The largest gains over libdeflate 1.26, one core, on the BGZF blocks of public genomics files
(an Intel Xeon Gold 5418Y, and an Apple M4 Pro where noted):

| | Data | zflate / libdeflate |
|---|---|---:|
| Decompression | bgzipped VCF | 1.7–1.9× |
| | long-read BAM, M4 Pro, two blocks at a time | up to 1.5× |
| | long-read BAM | up to 1.26× |
| | BCF | up to 1.24× |
| Compression, level 1 | BCF | 2.2–2.3× |
| Compression, level 5 | VCF | up to 1.8× |
| | BCF | 1.7× |

[`docs/BENCHMARKS.md`](docs/BENCHMARKS.md) has every result, per file and level, including where
zflate is slower; [`docs/PERFORMANCE-LOG.md`](docs/PERFORMANCE-LOG.md) records how they were
measured and every change tried on the way.

## Safety and testing

The hot loops use `unsafe` for unchecked loads and SIMD intrinsics, each covered by a
`SAFETY` comment that gives its proof. The checks:

- **Round trips** at every level and in every format through zflate, libdeflate and zlib.
- **Differential tests against libdeflate**: the same accept/reject decision on about 650,000
  damaged streams, and the same output on real blocks.
- **Guard pages** around inputs and outputs.
- **AddressSanitizer** and **Miri** runs.
- **Fuzz targets** for both directions (`fuzz/`, with cargo-fuzz).
- **Every block of eleven whole BAMs** (16.4 million blocks) decoded and checked against
  its CRC.

Tests that need real BGZF files read them from `ZFLATE_CORPUS`, a `:`-separated list of
BAM, BCF or bgzipped files; without it, they skip that part. The benchmarks (`cargo bench`)
need it too.

## Minimum Rust version

1.91.

## Licence

MIT ([`LICENSE`](LICENSE)). zflate is derived from libdeflate, whose MIT licence
is in [`NOTICE`](NOTICE).
