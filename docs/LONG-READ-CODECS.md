# Long-read codec comparison: x86, level 5

2026-11-07. Single-core **x86-64** testing, focused on zflate level 5, on two
aligned ONT samples and one unaligned HiFi sample. This is a separate run from
the layout-robust benchmarks in the [benchmark overview](BENCHMARKS.md).

Tested zflate commit: `6b61f30bfe787a87348d5d3eff66708f01d95540`.
The subsequent changes through `45d8d05293d6c87ecce65dd62a68afebdf7bb75f`
were reviewed: documentation, packaging, CI, benchmark corpus defaults and
equivalent lifetime elision in a CRC helper. No compression/decompression
algorithm changed. The results remain attributed to the actual tested commit.

## Level-5 compression

Seven timed rounds after one warmup, using **every block of the bounded samples**:
4,549 blocks and 294,998,192 uncompressed bytes. Speeds are uncompressed MiB/s.
ONT combines its two samples by summing their bytes and median elapsed times.
Size change is relative to libdeflate level 5; positive means larger output.

| Long-read workload | zflate L5 MiB/s | libdeflate L5 MiB/s | zflate speed advantage | zflate size change |
|---|---:|---:|---:|---:|
| ONT BAM | 146.4 | 122.7 | +19.3% | -0.13% |
| ONT CRAM 3.0 | 110.0 | 83.5 | +31.7% | -1.00% |
| ONT CRAM 3.1 | 128.1 | 86.4 | +48.2% | +1.91% |
| HiFi BAM | 173.0 | 146.7 | +17.9% | +0.24% |
| HiFi CRAM 3.0 | 148.2 | 97.6 | +51.9% | +8.83% |
| HiFi CRAM 3.1 | 148.2 | 97.7 | +51.7% | +8.87% |

**CRAM rows measure only gzip-coded data blocks**, excluding headers. They are
not whole-file CRAM sizes or end-to-end CRAM reader/writer timings. CRAM 3.0 and
3.1 were generated from the same selected BAM records with samtools 1.24.

## Speed at comparable size

Numeric compression levels do not mean equivalent work or output size across
libraries. The survey therefore compares zflate L5 with competitor settings at
0%, 1% and 2% output-size allowances. It includes levels 2, 3, 4 and 6 of
libdeflate, zlib-ng, zlib-rs and zenflate alongside the initial survey settings.

Within a **1% larger-output allowance**, zflate L5 was the fastest DEFLATE setting
tested on every workload but HiFi CRAM's gzip blocks. With **2%**, libdeflate L3
and ISA-L L3 also overtake it in places. On every block of the bounded samples
(the `full-confirm` phase of [summary.csv](results/x86-long-read-codecs-2026-11-07/summary.csv)),
these DEFLATE settings were faster than zflate L5:

| Workload | Faster DEFLATE setting | Speed relative to zflate L5 | Output relative to zflate L5 | Within |
|---|---|---:|---:|---:|
| HiFi CRAM 3.0 gzip blocks | libdeflate L3 | 1.110× | +0.36% | 1% |
| HiFi CRAM 3.1 gzip blocks | libdeflate L3 | 1.102× | +0.36% | 1% |
| ONT CRAM 3.0 gzip blocks | ISA-L L3 | 3.000× | +1.71% | 2% |
| HiFi BAM | libdeflate L3 | 1.045× | +1.29% | 2% |
| ONT CRAM 3.1 gzip blocks | libdeflate L3 | 1.034× | +1.41% | 2% |
| ONT BAM | libdeflate L3 | 1.023× | +1.49% | 2% |

The [matched-size table](results/x86-long-read-codecs-2026-11-07/matched-size.csv)
lists the best setting at each allowance from the up-to-8-MiB survey and sweep
phases, which cover more settings; its figures differ slightly from the
full-sample ones above (libdeflate L3 on HiFi CRAM 3.0: 1.096× at +0.55%), and
it agrees on where zflate L5 leads. Cross-phase comparisons use each phase's own
zflate L5 baseline. Differences of a few percent should be treated as near-ties,
not proof of a statistically significant lead.

These are results for these long-read workloads, not a general ranking of
codecs.

## Decompression

Twelve DEFLATE decoder implementations received identical payloads from four
reference encoders: libdeflate L5, ISA-L L3, zlib-ng L5 and zflate L5. ISA-L has
levels 0–3; decompression itself has no level setting. Each of the nine
sample/format cases was checked against all four references.

zflate had the highest median throughput in **31 of 36** comparisons. The five
other highest medians were on HiFi CRAM streams, with leads of approximately
1–5% for libdeflate or ISA-L. This is a result for the tested x86 corpus, not a
universal decoder ranking. Own-stream decoding results are recorded separately
and are not equal-input comparisons.

fdeflate receives the same raw DEFLATE payload wrapped in its required zlib
header and Adler32 trailer, so its framing/checksum work remains included.
The other DEFLATE decoder comparisons use raw payloads.

## Other formats and genomic codecs

The survey also includes C and Rust implementations of Zstandard, LZ4, Snappy,
Brotli, bzip2 and XZ; HTSCodecs and Noodles implementations of rANS, arithmetic,
FQZComp and read-name tokenization. These encodings are not interchangeable with
BAM's DEFLATE payloads.

Examples from the up-to-8-MiB survey, relative to zflate L5:

| Payload | Alternative | Speed | Output-size change |
|---|---|---:|---:|
| ONT CRAM 3.0 gzip-block contents | Zstandard C L3 | 3.15× | −0.12% |
| HiFi CRAM 3.0 gzip-block contents | Zstandard C L1 | 7.47× | −46.70% |
| HiFi CRAM 3.1 gzip-block contents | Zstandard C L1 | 7.83× | −47.74% |

Specialized FQZComp/read-name codecs produced smaller output while zflate L5 was
faster on those stream tests. For example, the smallest tested specialized
quality output was 24.5% smaller on ONT and 29.0% smaller on HiFi than zflate L5.
fdeflate often expanded these genomic bytes; that throughput is not an
equivalent-compression win and can require different BGZF chunking.

## Corpus and method

- Only records with stored sequence at least 1,000 bases were selected: up to
  2,000 records or 100 MiB of SAM per sample. These are bounded real-data samples.
- ONT HG003: 1,293 records, median 14,031 bases; ONT HG002: 1,414 records, median
  20,351 bases; unaligned HiFi HG002: 2,000 records, median 16,914.5 bases.
- CRAM 3.0/3.1 used `store_md=1,store_nm=1`. Record round trips checked all fields
  and auxiliary tags with tag order normalized and `decode_md=0`; undefined
  unmapped MAPQ was normalized for the HiFi validation comparison only. The
  benchmark BAM records were not changed by that comparison normalization.
- Gzip blocks account for about 41–43% of ONT CRAM 3.0 file bytes, 1.45–3.11% of
  ONT CRAM 3.1, and 60–62% of these unaligned HiFi CRAMs. These are **file-size
  fractions, not runtime fractions**.
- The initial survey, additional Rust implementations and nearby-setting sweep
  retain original independent block boundaries. Selection is deterministic,
  up to 8 MiB per case and including its largest block. Full confirmation keeps
  all source blocks. Read qualities use 100-record chunks up to 8 MiB; names
  retain NUL terminators and include all selected records.
- One warmup and seven timed rounds; rotating codec order; medians. Contexts and
  output buffers are reused where APIs permit. File I/O, initial setup and
  correctness checks are excluded; resets, required allocations and adapter
  copies are included. FQZComp C receives true read lengths and BAM flags;
  Noodles exposes a lengths-only FQZComp API.
- A process is pinned to one core. Dynamic clocks were not locked. Up to three
  attempts are allowed for interference: process utilization below 97% or
  sibling utilization above 10%. Attempts are selected by interference, not
  by the fastest codec timing. Per-sample median absolute deviations are retained.
- Rust 1.98.1, release, fat LTO, one codegen unit, default target with library
  runtime dispatch; optimized static native builds. This is one build per
  harness phase, not a layout-robust multi-build experiment.
- No accelerator, multicore, peak-memory or end-to-end pipeline measurements.
  No long-read BCF/VCF corpus was added.

## Correctness and coverage

All retained configurations passed exact-byte round trips, with independent
C/Rust cross-decoding where available. Original gzip members were checked by
both zflate and libdeflate, including their checksums. No zflate correctness
failure occurred in this comparison.

Some Noodles rANS-Nx16 order-1 configurations and arithmetic configurations
failed validation and were excluded from rankings (45 configuration/case
exclusions). Failures included encoding/decoding errors and byte mismatches;
some mismatches were also observed with the independent C decoder. This records
the harness outcome, not a completed upstream diagnosis. Noodles 0.102.0 required
visibility-only changes to expose private codec functions; no algorithms changed.

`deflate-rs` is encode-only. Its output was validated with libdeflate; no
libdeflate decoding speed is attributed to `deflate-rs`. Partial-coverage codec
results must not be compared as if they covered both ONT samples: use the
`samples` and `expected_samples` columns.

## Versions

| Family | Implementations tested |
|---|---|
| DEFLATE | zflate at the commit above; libdeflate 1.26 via libdeflater 1.26.1; zlib 1.4; zlib-ng 2.3.3; ISA-L 2.32.1; zlib-rs 0.6.8; miniz_oxide 0.9.1; zenflate 0.4.1; fdeflate 0.3.7; zune-inflate 0.2.54; libflate 2.3.2; deflate 1.0.0; yazi 0.2.2 |
| Zstandard | C 1.5.7 via zstd 0.13.3; ruzstd 0.9.1, implemented Fastest preset |
| LZ4 | C 1.10.0 via lz4-sys 1.11.1; lz4_flex 0.14.0 |
| Snappy | C++ 1.1.10 through its C API; snap 1.1.2 |
| Brotli | C 1.1.0; Rust brotli 9.1.1 |
| bzip2 | C 1.0.8; libbz2-rs-sys 0.2.6 with its Rust backend |
| XZ | liblzma C 5.4.5; lzma-rust2 0.25.0 |
| Genomic | HTSCodecs 1.6.7; Noodles CRAM 0.102.0 |

`zenflate-effort` is the same zenflate implementation using its native effort
presets. `deflate-rs` labels 1 and 6 mean Fast and Default. A label of 0 sometimes
means a library's default preset. These labels do not imply comparable effort.

## Recorded results

17,493 selected timing measurements and 2,499 per-sample setting summaries:

- [Level-5 confirmation](results/x86-long-read-codecs-2026-11-07/level5-comparison.csv).
- [All phases, pooled summaries](results/x86-long-read-codecs-2026-11-07/summary.csv).
- [Matched-size comparisons](results/x86-long-read-codecs-2026-11-07/matched-size.csv).
- [Per-sample medians and dispersion, gzip CSV](results/x86-long-read-codecs-2026-11-07/samples.csv.gz).
- [Individual timing measurements, gzip CSV](results/x86-long-read-codecs-2026-11-07/measurements.csv.gz).
- [Validation exclusions](results/x86-long-read-codecs-2026-11-07/validation-exclusions.csv).
- [Architecture, phase settings, input sizes and hashes](results/x86-long-read-codecs-2026-11-07/metadata.json).
- [SHA-256 checksums](results/x86-long-read-codecs-2026-11-07/SHA256SUMS.json).

Each phase is labeled; do not average the full-sample and up-to-8-MiB speeds
together. Raw measurements retain round and attempt identifiers; local core IDs
are omitted. Input hashes identify the sampled data; large genomic inputs are
not included. This folder contains results and methodology, not a self-contained
copy of the external comparison harness or its native dependency builds.
