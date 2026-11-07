# zflate: a pure-Rust DEFLATE decoder against libdeflate — performance log

*Written 2026-11-06 and updated with the first review's suggestions (§5.8), the long-read
work (§5.9), a layout finding (§5.10), an attempt at 150% (§5.11), the second review's
experiments (§5.12) and Apple silicon (§5.13, §5.14). It records what failed as well as what
worked, with the numbers. Commit hashes in these documents refer to zflate's development history before its public release; the public repository starts from that code as one commit.*

## 1. Summary

**Goal.** pasta (BAM/BCF/CRAM/VCF tooling) calls libdeflate through `libdeflate-sys` for BGZF
block compression and decompression, CRAM gzip blocks, and CRC32, because the Rust
alternatives were slower. The goal was a Rust-native codec that matches or beats libdeflate,
later raised to 120% of libdeflate's decompression speed, then to beating it as far as
possible on long-read BAMs (ONT, PacBio HiFi).

**Where it settled.**

- **CRC32**: own implementation (VPCLMULQDQ/AVX-512, PCLMULQDQ, slice-by-8, runtime
  dispatch). Equal to libdeflate at 4 KiB–64 KiB (the BGZF range), ahead at ≤16 bytes,
  behind by 5–15% at 64 B–2 KiB, 13% behind at 1 MiB.
- **Decompression (raw DEFLATE, zlib, gzip)**: own decoder, faster than libdeflate on every
  one of 22 BGZF files tested (§3.1), tuning and held-out:
  - **long-read BAMs (ONT, HiFi): 106–126%, 112.6% over the held-out set** (from 95–99%
    before the long-read work);
  - short-read BAMs: 105–109%;
  - BCF: 106–124%, 115.5% over the held-out set;
  - bgzipped VCF: 170–191%.
  With the block CRC check included, as pasta uses it: 105–123% on BAM, 108–120% on BCF.
  Byte-for-byte identical accept/reject behaviour to libdeflate on all tested input,
  including ~650,000 damaged streams.
- **What got long reads there**: entries that decode two literals in one lookup, built
  during the table fill, in a 12-bit table, only for streams that compress 8:1 or less
  (§5.9).
- **120% overall was not reached**: it was on cohort BCFs, the Arabidopsis HiFi BAMs and VCF,
  not on ONT, human HiFi, Illumina or the remaining BCFs. §6 says why.
- **Layout matters**: where the linker places the hot functions moves results by up to
  4–6% (§5.10); building pasta with `-C llvm-args=-align-all-functions=6` pins it.
- **150% was not reached, and looks out of reach on this CPU** (§5.11): a single stream
  is close to its floor (lookup latency, mispredicted literal/match branches, table
  building), and a decoder of several BGZF blocks at once without those branches lost on
  every file. Nothing from that attempt is in the crate.
- **The second review's experiments all lost** (§5.12): literal + length entries (HiFi
  108% → 98%), fused length + offset entries built during the fill (BCF 115% → 86%), and PGO
  (no gain over the function-alignment flag). ISA-L, as an independent reference, is slower
  than libdeflate on BAM and BCF except ONT (52–57% on BCF), faster on VCF (141–145%), and
  slower than zflate on every file. The code is unchanged.
- **On whole real files** (§3.3): every block of eleven whole BAMs (16.4 million blocks,
  900 GB inflated, including GIAB Illumina, HiFi and ONT whole genomes) and every gzip
  block of two CRAMs decoded correctly; timed per core with CRC32 against libdeflate on the
  same blocks: Illumina 103–104%, ONT 106%, human HiFi 113%, other HiFi 110–121%, CRAM gzip
  blocks 119–124%.
- **Validation** (§3.5): every `unsafe` site reviewed; the whole suite also run on the
  portable (non-BMI2, non-SIMD) code, under AddressSanitizer, and against guard pages
  (37,962 real blocks, 1,080,216 damaged streams); the decoder also run under Miri. That
  found one difference from libdeflate, now fixed: an error's kind on truncated input with
  exactly one byte of output room. No memory error or undefined behaviour was found.
- **Apple silicon** (§5.13, M4 Pro): one block at a time 115–121% on long-read BAMs;
  `deflate_decompress_batch`, two blocks at once, 142–150% with the CRC check. On Arm
  the CRC32 is 7–17% ahead of libdeflate from 2 KiB up.
- **Since done**: the compressor, all 13 levels (`COMPRESSOR-PLAN.md`), and pasta's switch
  from libdeflate to zflate for decompression, compression and CRC32, with the batch
  decoder in its multithreaded BGZF reader on Arm.

## 2. Setup and how things were measured

**Machine.** Intel Xeon Gold 5418Y (Sapphire Rapids), virtual machine, 24 vCPUs, AVX-512,
VPCLMULQDQ, BMI2. The VM exposes **no hardware performance counters** (`perf` reports
`cycles:u` unsupported), so profiling used timer sampling (`perf record -e cpu-clock`) and
experiments; no cycle or branch-miss counts were available. Rust 1.99.1 stable; no nightly.

**Corpus.** pasta's release-benchmark inputs, split into their BGZF blocks (≤64 KiB each):
`x.bam` (4,734 blocks), `c3.bcf` (6,760 blocks) for tuning, and for validation
`short.name.bam`, `fixed.bam`, a second `x.bam`, `hi.bcf`, `sa.bcf`, `lo.bcf`, `c4.bcf`, and
`fragmented-reference.fa.gz`. Throughput is uncompressed bytes per second, one core.

From the review round on (§5.8), a wider corpus: 48 MB windows of whole BGZF blocks from the
middle of each file (or at a set fraction of it), in tuning and held-out sets:

| Set | Files |
|---|---|
| Long reads, tuning | GIAB HG001 ONT (dorado), window at 30%; GIAB HG002 PacBio HiFi chr20; Arabidopsis HiFi (minimap2) |
| Long reads, held out | the same ONT run at 70%; GIAB HG002 Revio whole genome at 25%; fly HiFi (bonefish); fly dense reads; Arabidopsis HiFi (bonefish); fly reads vs. targets |
| BAM, held out | ONT HG001, Arabidopsis HiFi (minimap2), fly HiFi, pasta-written Illumina, short.name.bam |
| BCF, held out | cohort la.bcf and sorted.bcf; sa, lo, c4 |
| VCF.gz | 1000 Genomes high-coverage chr22 panel; cohort x.vcf.gz |
| Mix (regression check) | x.bam, Illumina chr20, c3.bcf, cohort x.bcf, 1000 Genomes chr22 |

Only public reference samples (GIAB) and non-human data. The BAM held-out set predates the long-read work, which
then tuned on two of its files (the ONT run, at another window, and the Arabidopsis
minimap2 BAM), so it is no longer fully held out; the long-read held-out set is.

The machine is shared (assembly and alignment jobs ran during the work), so every run waits
for a core that is under 10% busy and is pinned to it.

**Methods, in the order they were adopted.**

1. Criterion benchmarks (`benches/codec.rs`, `benches/crc32.rs`) — medians, one binary.
2. A profiling example (`examples/profile.rs`) run best-of-5, pinned to one core with
   `taskset`; zflate and libdeflate run back to back in the same script.
3. A **layout-robust** comparison (late in the work): each candidate built three times with
   `-C llvm-args=-align-loops={16,32,64}`, best-of-3 per build, ratio averaged across the
   builds. Adopted once it was clear that single builds vary by more than the effects
   being measured (next section).

**Measurement caveats that matter for interpreting every number here.**

- **Code layout moves zflate's results by ±1.5% on BAM and up to ±4% on BCF** between builds
  of identical source (BCF ranged 105.5%–113.9% of libdeflate across the three loop
  alignments). libdeflate (a C library linked into the same binary) stayed within ±0.3%.
  Single-build differences under ~2% are therefore not trustworthy.
- **Machine load was checked and is not the cause**: load average 8–10 of 24 cores
  (other jobs used ~2.5 cores); pinned and unpinned runs agreed, and libdeflate's numbers
  were stable all day.
- Criterion and the profiling example are different binaries with different layouts, so
  they disagree by a point or two; both are reported where it matters.

## 3. Final state

### 3.1 Decompression, every corpus file

Final code, profiling example, best of 3 per build, each run pinned to an idle core. The
zflate / libdeflate column is the mean over three builds with loop alignment 16, 32 and 64
bytes. The rest of the layout is wherever the linker put it. The last column repeats that
with every function aligned to 64 bytes (§5.10). The GiB/s columns are means over the three
builds. Files over 48 MB are measured on a 48 MB window of whole blocks.

| File | Set | libdeflate GiB/s | zflate GiB/s | zflate / libdeflate | with aligned functions |
|---|---|---|---|---|---|
| GIAB HG001 ONT, window at 30% | long read, tuning | 0.372 | 0.396 | 106.5% | 106.7% |
| GIAB HG002 HiFi chr20 | long read, tuning | 0.757 | 0.811 | 107.2% | 107.3% |
| Arabidopsis HiFi (minimap2) | long read, tuning | 0.647 | 0.809 | 125.1% | 125.3% |
| GIAB HG001 ONT, window at 70% | long read, held out | 0.374 | 0.400 | 106.8% | 107.0% |
| GIAB HG002 Revio HiFi, window at 25% | long read, held out | 0.890 | 0.960 | 107.9% | 108.0% |
| fly HiFi (bonefish) | long read, held out | 0.402 | 0.456 | 113.5% | 113.7% |
| fly dense reads | long read, held out | 0.403 | 0.455 | 112.8% | 113.0% |
| Arabidopsis HiFi (bonefish) | long read, held out | 0.648 | 0.814 | 125.6% | 125.7% |
| fly reads vs. targets | long read, held out | 0.370 | 0.419 | 113.2% | 113.6% |
| x.bam (pasta benchmark) | BAM, tuning | 2.862 | 3.010 | 105.2% | 105.7% |
| GIAB Illumina chr20 | BAM, tuning | 0.802 | 0.864 | 107.7% | 107.8% |
| GIAB HG001 ONT, window at 50% | BAM, held out | 0.400 | 0.426 | 106.5% | 106.8% |
| Illumina, written by pasta | BAM, held out | 0.786 | 0.843 | 107.3% | 107.5% |
| short.name.bam | BAM, held out | 3.204 | 3.480 | 108.6% | 108.8% |
| c3.bcf | BCF, tuning | 3.195 | 3.595 | 112.5% | 113.0% |
| cohort x.bcf | BCF, tuning | 2.585 | 3.193 | 123.5% | 124.1% |
| hi.bcf | BCF, tuning | 3.175 | 3.564 | 112.3% | 112.9% |
| cohort la.bcf | BCF, held out | 2.590 | 3.189 | 123.1% | 123.3% |
| cohort sorted.bcf | BCF, held out | 2.595 | 3.200 | 123.3% | 123.7% |
| sa.bcf | BCF, held out | 2.896 | 3.081 | 106.4% | 106.7% |
| lo.bcf | BCF, held out | 3.116 | 3.500 | 112.3% | 112.4% |
| c4.bcf | BCF, held out | 3.156 | 3.539 | 112.2% | 112.6% |
| 1000 Genomes chr22 VCF | VCF | 1.636 | 2.780 | 169.9% | 170.2% |
| cohort x.vcf.gz | VCF | 1.344 | 2.571 | 191.4% | 191.6% |

Over each set (time-weighted): long reads 111.3% tuning, **112.6% held out**; BAM 107.1%
tuning, 112.0% held out; BCF 116.4% tuning, **115.5% held out**; VCF 181.1%. Before the
long-read work the same held-out sets measured 97.2% (long reads), 97.1% (BAM) and 113.9%
(BCF).

Two older corpus files, one 64-byte-aligned build: bgzipped FASTA
(`fragmented-reference.fa.gz`) **120.8%** (was 100.4%; its 2–3-bit nucleotide codes pair
easily), `fixed.bam` 107.9%.

Other Rust decoders (criterion, earlier code): zune-inflate 1.98 / 2.30 GiB/s, zlib-rs
1.56 / 1.68 GiB/s on x.bam / c3.bcf, against libdeflate's 2.86 / 3.25.

### 3.2 pasta's actual pipeline: decompress, then verify the block CRC

One plain build, best of 3:

| File | libdeflate + its CRC | zflate + its CRC | |
|---|---|---|---|
| GIAB HG001 ONT | 0.369 GiB/s | 0.391 | 106.0% |
| GIAB HG002 HiFi chr20 | 0.747 | 0.791 | 105.9% |
| Arabidopsis HiFi | 0.640 | 0.788 | 123.1% |
| x.bam | 2.709 | 2.845 | 105.0% |
| Illumina chr20 | 0.791 | 0.843 | 106.6% |
| c3.bcf | 2.982 | 3.228 | 108.2% |
| cohort x.bcf | 2.459 | 2.954 | 120.1% |
| 1000 Genomes chr22 VCF | 1.571 | 2.601 | 165.6% |

### 3.3 Whole real files: every block checked, and timed

`examples/validate.rs` reads a whole file, checks every BGZF block's output against the CRC32
and length in its trailer, or every gzip block of a CRAM 3 file against libdeflate's output
(after checking each CRAM block's own CRC32, which confirms the parse), on several threads.
With `-time`, each worker then decodes every 16 MB chunk of blocks once more with each
decoder, on the same core, alternating which goes first, with the CRC32 check as pasta
does it; the figure is per core, four workers at once on a shared machine, so the GiB/s
are below the single-threaded ones of §3.1. Public data only: GIAB samples and non-human
genomes.

| File | Blocks | Compressed → inflated | zflate / libdeflate, with CRC32 |
|---|---|---|---|
| Illumina chr20 (as in §3.1) / the same written by pasta | 112,123 / 111,414 | 1.09 → 7.30 GB / 1.11 → 7.30 GB | 103.8% / 103.1% |
| Illumina HG003 NovaSeq 35×, whole genome | 4,852,927 | 46.80 → 315.83 GB | 104.4% |
| PacBio HiFi HG002 chr20 | 120,041 | 1.52 → 6.00 GB | 112.6% |
| PacBio Revio HG002, whole genome | 5,351,415 | 67.27 → 267.52 GB | 113.1% |
| ONT HG001, whole genome | 5,097,827 | 115.57 → 254.64 GB | 105.7% |
| Arabidopsis HiFi (minimap2 / bonefish) | 146,559 / 149,839 | 1.63 → 7.61 / 1.68 → 7.82 GB | 121.2% / 120.3% |
| fly HiFi, dense reads, reads vs. targets | 208,230 / 208,307 / 61,401 | 5.08 → 11.42 / 5.07 → 11.42 / 1.61 → 3.35 GB | 110.6% / 109.8% / 111.0% |
| CRAM 3.1, Illumina chr20: its gzip blocks | 19,246 of 63,910 | 0.06 → 0.30 GB | 119.2% |
| CRAM 3.0, ONT HG003 R10.4.1, whole genome: its gzip blocks | 663,644 of 1,201,758 | 2.47 → 9.74 GB | 123.5% |

**Every block was correct**: 16,420,083 BGZF blocks (248 GB compressed, 900 GB inflated)
and 682,890 CRAM gzip blocks, no failure. The whole files give lower ratios than the 48 MB
windows for Illumina (103–104% against 107%, as the whole files decode faster on average,
so table building weighs less) and higher ones for HiFi (113% against 107–108%).

### 3.4 CRC32 alone (GiB/s)

| Bytes | libdeflate | zflate | crc32fast | zune's slice-by-8 (start) |
|---|---|---|---|---|
| 16 | 3.7 | 3.8–4.0 | 2.7 | 3.1 |
| 64 | 12.2 | 10.4 | 7.6 | 2.3 |
| 128 | 22.3 | 20.7 | — | — |
| 256 | 34.1 | 33.1 | 20.1 | 1.6 |
| 1 KiB | 58.4 | 55.8 | 32.7 | 1.4 |
| 4 KiB | 71.2 | 70.2 | 55.1 | 1.3 |
| 64 KiB | 66.1 | 66.1 | 60.4 | 1.3 |
| 1 MiB | 76.6 | 66.7 | 60.7 | 1.3 |

### 3.5 Correctness

- Differential tests against libdeflate (`tests/differential.rs`): round trips at every
  level 0–12 in raw DEFLATE, zlib and gzip; exact, short, and oversized output buffers;
  damaged streams (bit flips, truncations, byte replacements) must be accepted exactly when
  libdeflate accepts them, with identical output — run up to 25× (~650,000 streams); every
  real BGZF block with its CRC; every output size around the fastloop's room limit and every
  truncation of the last 64 input bytes; every CRC32 implementation at every length to
  2 KiB, offsets 0–63, and boundary lengths to 1 MiB.
- With literal pairs: the same tests (the samples include streams under the 8:1 gate, so
  pairs are built and decoded in both loops), plus 9,950 real long-read and Illumina BGZF
  blocks from ten files, each checked against its block CRC.
- A unit test runs the non-BMI2 decoder build directly, since tests on this CPU otherwise only
  exercise the BMI2 build.
- Bugs found along the way: zune-inflate (the starting point) accepted zlib streams with a
  preset dictionary and code-length repeats that overrun the table, both of which libdeflate
  rejects; own bugs caught by the tests were a too-strong debug assertion at end of stream and
  an underflow in the bit-reverse lookup at the subtable hand-off.
- **Guard pages** (`tests/guard_pages.rs`): each input and output buffer sits against an
  inaccessible page, once at its end and once at its start, so that any access outside it
  faults. 37,962 real blocks from 12 files (long reads, Illumina, BCF, VCF), and 1,080,216
  damaged streams (three inputs, six lengths, four levels, three formats, 5,000 variants
  each) at four output sizes (exact, one byte short, half, 300 bytes over), each compared
  with libdeflate: no fault and identical results, on the default and the portable build.
  A control that reads one byte past a buffer faults.
- **Portable code**: the testing-only feature `force-portable` skips the BMI2 and SIMD
  builds; the whole suite passes on it too, with 5,000 damage trials per stream and the
  real blocks.
- **Found and fixed**: with literal pairs, a damaged, truncated stream that leaves exactly
  one byte of output room returned `InsufficientSpace` where libdeflate returns `BadData`:
  libdeflate writes the first literal, then refills before the second, and that refill
  fails on input that has run out. Both reject the stream; only the kind differed. A stream
  with pairs that fails for room is now decoded again without them, which gives
  libdeflate's answer; the hot loop is unchanged instruction for instruction, and the cost
  measured nil (c3.bcf 115.2% against 115.1%, HiFi 108.2% against 108.2%).
- **`unsafe` review**, every site:

  | Site | Invariant |
  |---|---|
  | BMI2 and SIMD dispatch (decoder, code-length reader, CRC32) | Called only once the CPU features are detected |
  | Fastloop refill (8-byte reads) and pointers | Runs only while 25 input and 301 output bytes remain, which bound one iteration's reads and writes |
  | Unchecked main-table lookups | The mask has at most 12 bits; the table has at least 4,096 entries (checked at compile time), all rewritten for each table |
  | Match copy | Offset checked against the bytes written; the room includes the copy's 39-byte overshoot |
  | Literal stores | Two bytes, within the fastloop's room |
  | CRC32 loads | Fixed-size array references, or a mask that selects only the input's bytes |

- **AddressSanitizer** (nightly 2026-11-05, `-Zsanitizer=address -Zbuild-std`): the whole
  suite, on the default and the portable build, with 400 damage trials per stream, 200
  guard-page trials, and 4 MB windows of eight corpus files: no errors. ASan's ODR check was
  turned off (`detect_odr_violation=0`): it stops on a duplicate global in libdeflate's own
  C code (`decompress_template.h`, included once per CPU variant), not in zflate.
- **Miri** (same nightly, Stacked Borrows): `tests/vectors.rs`, 434 streams compressed by
  zlib (`tests/data/make_vectors.py`: stored, fixed and dynamic blocks in all three formats,
  with and without literal pairs, 0–4,096 bytes, and two real 64 KiB blocks from public
  data), each into exact, longer and one-byte-short outputs, plus damaged variants: no
  undefined behaviour (27 minutes). Under Miri no CPU features are detected, so it runs the
  portable entry point: the same Rust as the BMI2 build, compiled without the feature.
- `cargo fmt`, `cargo clippy -D warnings` clean on both builds.

## 4. What zflate is

- `src/crc32.rs`, `src/crc32/{scalar,x86}.rs` — CRC32. Folding constants computed by `const fn`
  (checked against libdeflate's table); runtime choice of VPCLMULQDQ/AVX-512, PCLMULQDQ, or
  slice-by-8.
- `src/inflate.rs`, `src/inflate/{bits,table}.rs` — the decoder: libdeflate's design (one
  lookup per symbol with packed table entries, 64-bit bit buffer with branchless refill, a
  fastloop with no per-symbol bounds checks, a generic loop near buffer ends, BMI2 build
  chosen at runtime), plus the table-building changes in §5.5 and §5.7 and, for streams
  that compress 8:1 or less, literal-pair entries in a 12-bit main table (§5.9). That choice
  is made per call from the input and output sizes; callers need not do anything.
- `src/lib.rs` — libdeflate-shaped API: `Decompressor::{deflate,zlib,gzip}_decompress(data,
  out) -> Result<usize, DecompressionError>`, `crc32`, `adler32` (via `simd-adler32`).
  Decodes straight into the caller's buffer.
- Unsafe code: the fastloop (one region, raw pointers, guarded by its room check), the match
  copy, unchecked main-table lookups, three CRC loads, and the BMI2/feature dispatch calls;
  each has a `SAFETY` comment. Table building, header parsing, the generic loop and the
  zlib/gzip wrappers are safe code.

## 5. Chronological log

Percentages are zflate ÷ libdeflate decompression throughput; absolute numbers in GiB/s.
"Kept" / "reverted" refers to whether the change stayed in the code.

### 5.1 CRC32

| Step | Result |
|---|---|
| Baseline: zune-inflate's slice-by-8 | ~1.3 GiB/s at all sizes ≥1 KiB (≈50× slower than libdeflate) |
| Ported libdeflate's folding design (8 accumulators, VPCLMULQDQ/AVX-512 and PCLMULQDQ, Barrett reduction, masked loads for 4–15 bytes) | First build slower than slice-by-8: a `const fn` called in ordinary code runs at run time, recomputing `x^N mod G` (~20,000 loop steps) per call. Binding results to `const` items fixed it: parity at 4–64 KiB |
| Small sizes: removed a hidden `memcpy` (from `as_chunks` + `zip` loads of N lanes) | 64 B 8.3 → 9.8; 1 KiB 52.5 → 55.8 |
| Bounded the tail loop (const max lanes) | small gain at 16 B only |
| Straight-line stage cascade (no `Option<[__m512i; N]>` merges, which went through the stack) | 128 B 17.6 → 20.7; 256 B 30.7 → 33.1; no stack accesses left in the function |
| Remaining 64 B gap | ~0.9 ns/call with near-identical instruction streams; undiagnosable without counters |
| Not done | libdeflate's 64-byte pointer alignment for inputs >64 KiB (the 1 MiB gap) |

### 5.2 Decoder baseline and match copy

The vendored zune-inflate decoder was replaced by an own, deliberately plain port of
libdeflate's decoder (same table format, same accept/reject rules, all safe code).

| Step | x.bam | c3.bcf | |
|---|---|---|---|
| libdeflate | 2.86 | 3.25 | |
| zune-inflate (stock crate) | 1.99 | 2.31 | |
| Plain own decoder (byte-wise copies, one refill per symbol) | 1.57 | 0.66 | |
| 8-byte word copies, safe slices | 1.84 | 1.90 | kept |
| One range check, then raw-pointer copies | 2.07 | 2.17 | kept |
| libdeflate's unrolled copy (5 words; 4 for runs; 2 for short offsets) | 2.12 | 2.63 | kept |
| 16-byte words for offset ≥16 | 2.00 | 2.63 | reverted |

### 5.3 Fastloop

| Step | x.bam | c3.bcf | |
|---|---|---|---|
| A: split fastloop / generic loop, skip end and room checks | 2.11 | 2.63 | kept (no gain alone: checks were well predicted) |
| B: look up the next entry before refill and copy; up to 3 literals per refill | 2.48 | 3.03 | kept |
| BMI2 build of the decoder, chosen at runtime | 2.50 | 3.10 | kept |
| `#[inline(always)]` on the match copy (the BMI2 build had stopped inlining it) | 2.62 | 3.16 | kept |
| Unchecked fastloop accesses | 2.64 | 3.17 | reverted (<1%, not worth the `unsafe`) |

Profiling then showed table building at parity with libdeflate and the whole remaining gap
in the fastloop, with ~6% of its samples on register spills (libdeflate's hot loop has none).

| Step | x.bam | c3.bcf | |
|---|---|---|---|
| Local copy of the bit reader while reading code lengths | no change | | reverted |
| Track input as a shrinking slice (fewer live registers) | 2.58 | | reverted |
| Generic loop moved out of line | 90.1% | 95.6% | reverted |
| **libdeflate's "subtract the whole entry" bit-count trick** | 91.3 → 91.8% | 95.5 → 98.4% | kept |
| **Fastloop over raw pointers with its own local bit state** (4 live buffer values instead of 8, as libdeflate) | 97.0% | 99.1% | kept |

### 5.4 Data analysis (temporary instrumentation, removed afterwards)

| | x.bam | c3.bcf |
|---|---|---|
| BGZF blocks / DEFLATE blocks | 4,734 / 9,407 | 6,760 / 12,103 |
| Repeated code-length sets between blocks | 0% | ≤1.5% (64-entry cache) |
| Share of output bytes from literals | 2.1% | 1.2% |
| Average match length | 37.4 | 75.1 |
| Match lengths | 64% ≤31; 35% 41–160 | 32% 8–15; 23% ≥161 |
| Offsets | 76% ≥256; 12% 16–31 | 20% exactly 1 (runs) |
| Literal/length lookups needing a subtable | 0.04% | ~0% |

This ruled out caching tables and multi-literal lookups, and pointed at match decoding and
table building.

### 5.5 Getting past libdeflate

| Experiment | x.bam | c3.bcf | |
|---|---|---|---|
| E1: unchecked main-table lookups (now that registers were free) | 97.9% | 103.5% | kept |
| E2: reuse tables across blocks | — | — | ruled out (0% repeats) |
| E3: literal/length table 10 / 9 bits | 95.5 / 92.5% | 102.3 / 92.3% | reverted (11 stays) |
| E3: offset table 9 / 10 bits | 96.9 / 97.2% | 102.5 / 103.8% | reverted (8 stays); re-tested later, still worse |
| E4: 16-byte copies, offset ≥16, 32-byte rounds | 99.3% | 98.6% | reverted |
| E4: 3×16 / 4×16 rounds | 89.4 / 98.2% | 100.6 / 102.3% | reverted |
| E4: 32-byte then 64-byte rounds | 91.3% | 103.5% | reverted (confirmed real with forced loop alignment, likely store-forwarding failures) |
| E5a: count code lengths with per-length `filter().count()` | 82.7% | 85.1% | reverted (did not vectorize) |
| **E5a′: four interleaved histograms** | **100.3%** | **104.8%** | kept — first time ahead on both |
| **E5b: sort symbols by length in four independent quarters** | **101.7%** | **105.8%** | kept |
| **E5c: fill the main table from canonical codes + an 11-bit reverse-lookup table** (breaks the bit-reversed-increment dependency chain) | **104.6%** | **108.9%** | kept |
| Inline code-length reading into the BMI2 decoder | 103.3% | 109.6% | reverted (layout) |
| Local bit reader in code-length reading | 104.5% | 109.6% | kept |
| `zip`-based quarter loops (no index checks) | 104.9% | 110.2% | kept |
| BMI2 code-length reading via a generic parameter | 103.5% | 109.1% | reverted |
| E4d: 16-byte copies only for offset ≥32 (2× / 3× per round) | 106.2 / 106.0% | 105.1 / 106.9% | reverted (BAM gain, BCF loss) |
| E4d′: same for offset ≥64 | 106.0 / 105.8% | 103.6 / 105.3% | reverted |
| Require AVX2 in the fast build (no code change) | 104.3% | 108.2% | reverted |
| E4f: 32-byte AVX2 copies, offset ≥32 / ≥64 | 104.1 / 103.9% | 100.4 / 100.0% | reverted |

The key insight of this phase: libdeflate's table builder spends most of its time in
**store-to-load dependency chains** (`count[len]++` and `sorted[offset[len]++] = sym` when
neighbouring symbols share a length, and the bit-reversed codeword increment). Breaking those
chains took table building from ~13.7% to ~7.6% of BAM decode time.

### 5.6 Attempts at 120% (all reverted)

| Attempt | x.bam | c3.bcf | Why it lost |
|---|---|---|---|
| Branch-free code-length decoding (repeat count and extra bits in the precode entry, 16-byte splat stores) | 103.7% | 108.0% | cost was the lookup chain, not mispredicts |
| Unconditional 80-byte first copy round | 100.8% | 109.0% | extra stores on short matches |
| **Two-block API**: decode two independent BGZF blocks interleaved in one call | 92.3% | 93.9% | see below |
| Unconditional refill before the offset (remove a data-dependent branch) | 99.9% | 99.5% | puts the refill on the critical chain |
| Fixed 11-bit literal/length table (constant mask) | 104.6% | 105.9% | BCF's short codes pay for a bigger table |
| **Fused CRC**: fold each new 64-byte chunk of output into the CRC inside the fastloop | 2.64 vs 2.81 GiB/s (zflate+CRC) | 3.26 vs 3.23 | see below |
| Fused CRC, batched every 8 / 16 iterations | 2.62 / 2.64 | 2.94 / 3.07 | as above |

**Two-block API.** Hypothesis: one stream's fastloop is a latency-bound chain, so a second
independent stream would fill idle execution units. Implemented as an additive API
(`deflate_decompress_pair`), verified correct (differential tests including damaged inputs
and all real block pairs), but slower on all four files tried (87–99% of libdeflate, against
100–107% for single-block zflate). The profile showed both streams' loops running, so the
interleaving worked mechanically; the likely cause is that the two streams' unpredictable
branches share the branch predictor, plus register pressure on x86-64's 16 GPRs. Code saved
(not in the crate). It may behave differently on ARM64 (31 GPRs) or other cores.

**Fused CRC.** Decoding with the AVX-512 features enabled but the CRC computed afterwards ran
at 2.82 GiB/s, so the compiler features were harmless; folding inside the loop, per iteration
or in batches, cost ~6% of decode time — more than the ~5% a separate CRC pass costs. Code
saved (not in the crate).

**Copy-cost measurement.** Skipping match copies entirely (timing only, wrong output):
x.bam 2.99 → 3.16 GiB/s, c3.bcf 3.55 → 4.21. Copying is only ~5% of BAM's time and ~16% of
BCF's; the rest is Huffman decoding.

### 5.7 Exhaustive table-building round (layout-robust measurement)

Line-level profile of the starting point: table building 7.5% (counting ~1.1%, sort ~2.0%,
fill ~1.6%, prefix sums ~0.5%, libc `memmove` for table doubling ~0.7%), code-length reading
8.3% (precode lookup → shift → subtract chain, ~9–10 cycles per code length).

| Experiment | x.bam | c3.bcf | |
|---|---|---|---|
| Baseline | 103.8% (103.6 / 103.6 / 104.2) | 109.7% (108.2 / 110.1 / 110.9) | |
| **T1: refill once per three code lengths** (no per-symbol "refill?" branch) | **105.2%** (104.8 / 105.0 / 105.8) | 109.7% | kept |
| **T2: BMI2 build of the code-length reader, dispatched inside it** | **105.6%** | **111.2%** | kept |
| T3b: inline 16-byte stores for repeat runs instead of `memset` | 105.6% | 111.2% | reverted (neutral) |
| T3a: inline 64-byte-chunk table doubling | 104.8% | 110.9% | reverted (LLVM re-emitted `memcpy`) |
| T4: per-length bitmasks (SSE2 compares + popcount) replacing count + sort | 100.0% | 104.5% | reverted |
| T5: one consume per repeat code (extra bits in the precode entry) | 104.8% | 109.4% | reverted |
| T6: cache precode / offset tables when lengths repeat | — | — | ruled out: precode repeats 2.8% (BAM) / 4.6% (BCF), offset lengths 0% |

### 5.8 Review round: the reviewer's suggestions, and a wider corpus

The reviewer suggested multi-symbol entries: literal pairs, and length + distance fused into
one lookup. Eligibility measured on our files (temporary counters):

| | x.bam | c3.bcf | short.name.bam | sa.bcf |
|---|---|---|---|---|
| Literals, share of all symbols | 45.0% | 46.7% | 40.5% | 50.5% |
| Matches eligible for fusing | 84.2% | 48.3% | 82.8% | 52.3% |
| Literal pairs that fit one 11-bit lookup (share of symbols) | 2.3% | 0.7% | 0.07% | 0.4% |

A fused length + distance prototype (all differential tests passing) sped the decode loop
itself (profile share ÷ throughput, builder time excluded) by **~0.5% on x.bam** and **~12% on
c3.bcf**. By the reviewer's own stopping rule it was retired for BAM: BAM's loop is not limited
by the second lookup.

**The wider corpus changed the picture.** 48 MB windows from the middle of real files (GIAB
ONT and PacBio HiFi, Illumina, Arabidopsis and fly HiFi; more BCFs; 1000 Genomes VCF), split
into tuning and held-out sets, single build:

| Set | Files | zflate / libdeflate |
|---|---|---|
| BAM, tuning | x.bam 105.8%, Illumina chr20 **95.9%**, HiFi chr20 **99.5%** | 98.6% |
| BAM, held out | ONT HG001 **95.2%**, Arabidopsis HiFi **98.5%**, fly HiFi **97.3%**, pasta-written Illumina **96.1%**, short.name 109.4% | 97.1% |
| BCF, tuning | c3 110.3%, cohort x.bcf 123.0%, hi 110.4% | 114.8% |
| BCF, held out | la 122.3%, sorted 122.7%, sa 104.1%, lo 110.0%, c4 110.1% | 113.9% |
| VCF.gz | 1000 Genomes chr22 169%, cohort x.vcf.gz 191% | 181% |

x.bam, the file everything before had been tuned on, is not representative: it decodes at
2.9 GiB/s, real BAMs at 0.36–0.9 GiB/s, because quality strings and tags barely compress.
**On real BAMs zflate was 1–5% slower than libdeflate.** That led to the long-read work below.

### 5.9 Long-read BAMs (ONT, PacBio HiFi)

**Diagnosis.** Long-read BAMs compress only 2–5:1. 52–86% of their symbols are literals,
matches average 6–18 bytes, and decoding costs ~28 (HiFi) to ~32 (ONT) cycles per symbol,
against ~15 on x.bam. Both decoders spend ~95% of the time in the decode loop and table
building is only ~4%, so zflate's builder advantage does not show. The cost per symbol points
at mispredicted "literal or match?" branches, which no counters can confirm here.
An unconditional refill before the offset (removing one data-dependent branch) made it worse
(96.0% against 97.8%).

**The lever: literal pairs.** A main-table entry whose two codewords are both literals and
fit in the table's bits together decodes both with one lookup and one unconditional 2-byte
store (the output pointer advances by 1 or 2). This cuts the number of lookups and branches
on literal runs, which is exactly where long reads spend their time. Entry layout: `LITERAL |
PAIR`, first literal in bits 16–23, second in 8–15, both codewords' bits in 0–7, so the
existing whole-entry bit consumption works unchanged.

Steps, tuning set (three-layout means; Arabidopsis = minimap2-aligned Arabidopsis HiFi):

| Step | ONT | GIAB HG002 HiFi chr20 | Arabidopsis | Long-read total | Mix total¹ | |
|---|---|---|---|---|---|---|
| Baseline | 96.3% | 99.4% | 98.9% | 97.8% | — | |
| Pairs, built by scanning the table after the fill, gated in the hot function | 95.9% | 97.9% | 113.1% | 100.6% | 113.9% | |
| Pairs on every table (no gate) | 99.0% | 100.9% | 116.7% | 103.9% | — | lost on BCF |
| Gate decided in the dispatcher, outside the inlined decoder | 98.7% | 100.7% | 116.0% | 103.4% | 116.1% | kept |
| 12-bit main table where pairs are built (more pairs fit) | 104.3% | 100.4% | 112.4% | 105.4% | 117.8% | kept |
| 12 bits only where a code-length estimate predicts enough extra pairs | 103.8% | 99.8% | 115.3% | 105.7% | 117.1% | kept, later removed |
| + fused length + distance entries | 103.5% | **85.3%** | 105.9% | 99.1% | 113.1% | reverted |
| **Pairs written during the table fill**, before each doubling, instead of a scan | 106.2% | 105.5% | 123.6% | 110.4% | 118.2% | kept |
| Pairs on every table again, now that they are cheap | 106.2% | 105.8% | 123.9% | 110.5% | 116.1% | reverted (BCF/VCF lose) |
| **Always 12 bits where pairs are built** (estimate dropped) | **106.9%** | **106.5%** | **124.3%** | **111.1%** | **119.3%** | kept |
| 13-bit main table where pairs are built | 103.6% | 103.9% | 119.8% | 107.7% | 114.2% | reverted |
| **Final**: as above, with `build_tables` never inlined (§5.10) | **106.5%** | **107.2%** | **125.1%** | **111.3%** | — | kept |

¹ Mix: x.bam, Illumina chr20, c3.bcf, cohort x.bcf, 1000 Genomes chr22 VCF, to check that
nothing else regresses. Rows before the final one carry up to 4–6% of function-placement
noise that the three-build mean does not remove (§5.10); the close calls were re-checked
with a pinned layout.

What made it work:

- **Where to build pairs.** Pairs cost table-building time and pay back per symbol decoded,
  so they pay on poorly compressed streams (many symbols per table) and not on BCF/VCF
  (12–34:1, few symbols per table). The decoder knows the compressed and output sizes before
  it starts, so `inflate()` builds pairs when the ratio is at most 8:1. Deciding this inside
  the inlined hot function cost ~3 points on long reads by itself (layout and register
  pressure), so the flag is set once in the dispatcher.
- **How to build them.** The first builder scanned the finished table for literal entries and
  paired them: 31–37% of the time on long reads. The final one records each literal's
  bit-reversed codeword as the fill places it, and writes all pairs of a given total length
  just before the table doubles past that length, so the doubling copies them along like any
  other entry. Each pair slot is written once and nothing is scanned.
- **Table width.** At 11 bits, ONT's 6–7-bit quality codes rarely pair (~3% of symbols); at
  12 bits they do. 12 bits doubles the fill cost, which only pays where pairs are built; 13
  bits loses everywhere (fill cost and L1 pressure).
- **Fused matches lose on long reads** even on top of pairs: the larger entries' builder costs
  more than their loop saves when matches are short and few.

The generic (bounds-checked) loop decodes pairs too. Its subtable test became `entry &
(LITERAL | SUBTABLE_POINTER) == SUBTABLE_POINTER`, because a pair's second literal sits in the
flag bits. The fastloop's room check allows for its literal entries' 2-byte stores.

### 5.10 Function placement: the second source of layout noise

Removing the width estimate (dead code by then) made the build **4–6% slower** on long reads
and BCF alike: held-out BCF fell from 115.3% to 112.4%, the tuning set's long reads from
111.1% to 106.4%, while libdeflate's own speeds did not move. Two causes, both layout:

1. **Inlining.** With `build_tables` smaller, LLVM inlined one of its two call sites into the
   hot decoder, which then called `table::build` directly. As in §5.9, builder code inside
   the hot function costs the loop. `#[inline(never)]` on `build_tables` restored the long
   reads (ONT 100.3% → 105.9%).
2. **Function addresses.** The rest (c3.bcf 108.4% against 112.9%) followed where the
   linker placed the hot functions within 64-byte lines: the dead code shifted
   `inflate_bmi2` and `table::build` by 208 bytes. Built with every function aligned to 64
   bytes (`-C llvm-args=-align-all-functions=6`), the two versions measure the same, at the
   fast end: ONT 106.7% / 106.7%, c3.bcf 115.0% / 115.2%, x.bcf 125.8% / 125.7%.

The three-loop-alignment suite did not catch this, because all three builds place functions
the same way. From here on, decisions use an A/B of builds with 64-byte function alignment,
run alternately on one core, which repeats to about ±0.5%. Re-checked that way:

| Decision | ONT | HiFi chr20 | Arabidopsis | Illumina | x.bam | c3.bcf | x.bcf | VCF |
|---|---|---|---|---|---|---|---|---|
| Final: always 12 bits with pairs, ratio gate 8 | 107.0% | 108.3% | 126.0% | 107.8% | 106.2% | 115.5% | 125.6% | 171.1% |
| 12 bits only where the estimate predicts enough pairs | 106.5% | 106.9% | 124.8% | 106.5% | 106.4% | 115.4% | 125.5% | 171.1% |
| Pairs on every table (no gate) | 107.0% | 108.2% | 126.5% | 107.8% | 104.2% | 112.7% | 121.6% | 170.3% |
| 13-bit table with pairs | 105.9% | 106.7% | 122.6% | 102.1% | | | | |
| 2-word first copy round (then 5-word rounds) | 103.5% | 101.6% | 120.9% | 98.5% | 106.1% | 107.5% | | |
| 9-bit offset table | 102.2% | 103.6% | 121.1% | 103.4% | 103.9% | 110.6% | | |

So the earlier decisions hold. Stable Rust cannot align a function from inside a crate
(`#[repr(align)]` is rejected on functions, `rustc_align` is nightly-only), so pasta would get
the pinned layout only by building with `-C llvm-args=-align-all-functions=6`; otherwise its
speed depends on where zflate lands in its binary.

### 5.11 Aiming at 150%: what limits decoding, and a multi-block decoder

**Synthetic streams.** BGZF files written with a controlled token stream
(`synth.py`: dynamic Huffman, 7-bit literals, 20% matches of 4–12 bytes), so that costs
can be separated (GiB/s; cycles per token at the measured 2.78 GHz):

| Stream | libdeflate | zflate | Cycles per token |
|---|---|---|---|
| Literals only | 0.358 | 0.358 | 7.2 |
| Matches every fifth token, offsets 64–4096 | 0.639 | 0.645 | 8.9 |
| Same tokens, matches at random positions | 0.449 | 0.445 | 12.7 |
| Random, offsets 1–7 | 0.413 | 0.385 (93%) | |
| Periodic, offsets 1–7 | 0.568 | 0.518 (91%) | |

A literal/length lookup costs ~7 cycles (the load-shift-mask chain, the same in both
decoders); unpredictable literal-or-match decisions add ~3.8 cycles per token, about one
~17-cycle mispredict per match. zflate is behind on offsets under 8, but they are 3–5% of
matches in real long-read files.

**Real files** (token statistics from a Python inflate of 20–40 blocks each, `tokstats.py`):

| File | Bytes / token | Literals | Matches right after a match | Deflate blocks per BGZF block | zflate cycles / token |
|---|---|---|---|---|---|
| HiFi chr20 | 5.9 | 64% | 48% | 3.2 | ~19 |
| Illumina chr20 | 6.0 | 53% | 48% | 1.7 | ~18 |
| fly HiFi | 2.2 | 84% | 32% | | ~12 |
| Arabidopsis HiFi | 3.4 | 87% | 41% | 2.9 | ~11 |
| cohort x.bcf | 38 | 53% | 66% | 2.0 | ~25 (loop) |

Table building, by profile share: HiFi 15.6% (12-bit tables with pairs, three per BGZF
block, ~10,000 cycles each), c3.bcf 19.6%, x.bcf 20.1%, x.bam 15.8%, ONT 4.4%. So on HiFi
roughly 7 cycles per token are lookups, 6–7 mispredicts, 3 table building and 2 copies and
the rest; libdeflate takes 20.1. 150% would need 13.4.

**A multi-block decoder.** BGZF blocks are independent streams, so a
decoder working on two at once has two independent lookup chains. The first attempt at
this (§5.6) interleaved two copies of the usual loop and lost to the doubled mispredicts.
This one had no per-token branch: a first pass per stream always looks up an offset after
the literal/length entry and consumes its bits only for a length (a select), writes
literals to an array and each match's bits to a list; a second pass per block writes the
literal runs and matches out. Rare cases (subtables, end of block, the last input bytes)
branch out, and any failure re-decodes that block with the ordinary decoder. It decoded
every block of the test files correctly. Iterations (GiB/s; one block per call: HiFi
0.811, Arabidopsis 0.809):

| Version | HiFi chr20 | Arabidopsis | Why |
|---|---|---|---|
| First | 0.559 | 0.555 | LLVM compiled the select back into a branch; 66 instructions per step, with spills |
| `select_unpredictable`, a step budget instead of per-step end checks, lanes' tables side by side, 12-bit tables always | 0.708 | 0.681 | |
| Flag bits testable with one OR, `movbe` literal stores, `popcnt` advances, second pass in L1-sized chunks | 0.749 | 0.717 | first pass ~7.5 cycles per token, second pass ~20 per match |

The first pass reached what the design promised: two lanes at ~7.5 cycles per token,
against ~15 in the single decoder's loop. The second pass gave it all back: ~80
instructions per match (decode the match again, copy the literal run, copy the match)
with no decoding latency to hide behind. Skipping the match copies made the whole decoder
0.865 GiB/s, and copying from a buffer nothing had just written 0.796, so store-forwarding
stalls were only ~6%. On the other files it lost more (x.bam 1.91 against 3.01 GiB/s, c3.bcf
2.82 against 3.58). At ~30 instructions per symbol in the first pass, plus the second, the
issue width caps it near the single decoder; it was taken out of the crate.

**Other attempts, reverted:**

| Change | c3.bcf | x.bcf | x.bam | HiFi | ONT | Arabidopsis | VCF |
|---|---|---|---|---|---|---|---|
| Final code (aligned functions) | 115.1% | 125.9% | 106.4% | 108.2% | 106.7% | 125.8% | 171.1% |
| Precode pairs: two code lengths per lookup | 109.7% | 117.6% | 102.3% | 104.1% | 105.9% | 121.0% | 167.4% |
| Decoder split into header, block and block-loop functions (needed by the batch decoder) | 113.9% | | 105.7% | 104.9% | 102.4% | | |

Pairs in the precode table would cut code-length lookups by 29–38%, but that is worth only
~1.5% once the extra table pass is paid, and the version tried added a select to the
lookup chain. Fused length + offset entries for well-compressing streams (§7 question 1)
were not retried: their builder writes up to ~4,000 strided entries per table, more than
the ~12% loop gain on BCF saves; with a builder that writes each once before the table
doubles, as pairs are, the net would be a few percent.

### 5.12 Second review: a paired baseline, ISA-L, new entry kinds, PGO

The second review recommended one bounded round with explicit acceptance criteria (at least
2% whole-member throughput on the intended data, no reproducible regression elsewhere):
a paired, median-based baseline; ISA-L as an independent reference; fused length + offset
entries built at their natural depth for BCF; literal + length entries and arithmetic
extensions of literal pairs for BAM; and instrumentation-based PGO.

**Method.** One build with every function aligned to 64 bytes; each decoder run five times,
alternating order, each run pinned to an idle core; the figure is the median of the five
per-run ratios to libdeflate, with their range. ISA-L 2.31.1 (built here with nasm 2.16.03)
is called through `isal_inflate_stateless`, with `ISAL_GZIP_NO_HDR` for its own CRC.

**Baseline** (decode only; decode + CRC; ISA-L decode only; ISA-L decode + CRC):

| File | zflate | zflate + CRC | ISA-L | ISA-L + CRC |
|---|---|---|---|---|
| GIAB HG001 ONT, 30% / 70% | 106.7% / 106.9% | 106.8% / 106.7% | 101.3% / 101.6% | 100.8% / 101.1% |
| GIAB HG002 HiFi chr20 | 106.7% | 106.9% | 82.1% | 80.6% |
| GIAB HG002 Revio, 25% | 107.4% | 107.4% | 80.0% | 78.5% |
| Arabidopsis HiFi (minimap2 / bonefish) | 124.7% / 125.2% | 124.7% / 125.3% | 89.4% / 89.4% | 88.3% / 88.3% |
| fly HiFi, dense, reads vs. targets | 112.3–113.4% | 112.5–113.5% | 89.9–92.2% | 88.8–91.1% |
| Illumina chr20 / written by pasta | 107.4% / 106.8% | 107.2% / 106.9% | 87.5% / 87.6% | 87.3% / 87.0% |
| x.bam / short.name.bam | 106.3% / 108.9% | 106.0% / 108.2% | 72.8% / 71.5% | 71.8% / 70.7% |
| c3, hi, lo, c4 BCF | 111.8–112.5% | 111.2–111.9% | 56.4–57.0% | 57.4–57.7% |
| cohort x, la, sorted BCF | 123.1–123.4% | 122.0–122.6% | 51.9–52.2% | 52.7–52.8% |
| sa.bcf | 106.8% | 106.7% | 61.5% | 62.1% |
| 1000 Genomes chr22 VCF / cohort x.vcf.gz | 170.4% / 192.0% | 165.5% / 187.6% | 140.6% / 145.3% | 137.8% / 140.2% |

The ranges of the five runs were within about ±0.5% for most files. These agree with §3.1
and §3.2 to within a point. ISA-L gives no evidence of headroom on this CPU: it is
slower than zflate on every file, and slower than libdeflate on BAM and BCF except ONT
(it is ahead of libdeflate on ONT and VCF).

**Eligibility on the real windows** (`coverage.py`, 30 blocks per file: primary
literal/length lookups saved against today's greedy 12-bit pairs, each rule alone):

| File | Literal + length | Arithmetic pair extension | Both |
|---|---|---|---|
| Human ONT, 30% / 70% | 14.2% / 14.5% | 5.0% / 4.6% | 18.3% |
| Human HiFi chr20 / Revio | 15.0% / 15.4% | 0.3% / 0.2% | 15.3% / 15.5% |
| Illumina chr20 | 11.7% | 0.7% | 12.3% |
| fly HiFi | 8.8% | 8.1% | 16.0% |
| Arabidopsis HiFi | 6.8% | 0.0% | 6.8% |
| x.bam | 4.3% | 0.5% | 4.8% |

Fusable matches at an 11-bit root: x.bcf 73%, hi.bcf 64%, x.bam 90%.

**Literal + length entries** (a literal and the length after it, extra bits included, in one
main-table entry; built during the fill like pairs; decoded in the length path, which writes
the literal first). All tests passed, but:

| | ONT | HiFi chr20 | Revio | Arabidopsis | fly | Illumina |
|---|---|---|---|---|---|---|
| Baseline | 106.7% | 108.2% | 108.8% | 126.0% | 113.7% | 107.8% |
| Literal + length | 105.6% | 97.9% | 97.5% | 111.0% | 101.2% | 102.4% |

On HiFi the builder took 5.6% of the time, and the decode loop itself became 3.7% slower
despite 15% fewer lookups: the extra entry type costs the common path more than the saved
lookups.

**Fused length + offset entries, built at their natural depth** (each written once, before
the table doubles past its bits, for streams that compress better than 8:1; checked entry by
entry against a reference decode in a unit test):

| | c3 | x.bcf | sa | VCF chr22 | x.bam | short.name |
|---|---|---|---|---|---|---|
| Baseline | 115.1% | 125.6% | 108.2% | 170.8% | 106.1% | 109.0% |
| Fused | 86.4% | 95.0% | 83.3% | 147.4% | 90.6% | 100.5% |

On c3.bcf the decode loop became 17% slower (with 64–73% of matches fused, which kind an
entry is predicts badly), and the builder took 11.5% of the time. The ~12% loop gain of
§5.8 was measured on the pre-pairs code with profile shares; it does not reproduce.

**Arithmetic pair extensions** were not built: they would add work to every literal store,
and they save lookups only on ONT (5%) and fly (8%), not on HiFi; literal + length entries,
with three times the saving, already lost.

**PGO** (instrumentation-based, trained on the tuning files, measured on held-out files):

| File | A: plain | B: aligned functions | C: PGO | D: PGO + aligned |
|---|---|---|---|---|
| ONT, 70% | 106.1% | 107.2% | 107.0% | 107.2% |
| Revio | 106.5% | 108.2% | 107.4% | 107.6% |
| fly HiFi | 112.7% | 113.6% | 112.9% | 113.4% |
| Arabidopsis (bonefish) | 124.0% | 126.0% | 124.5% | 125.1% |
| Illumina, written by pasta | 106.2% | 107.4% | 107.5% | 107.8% |
| la.bcf / sorted.bcf | 121.0% / 121.2% | 124.1% / 124.3% | 120.6% / 120.8% | 121.5% / 121.6% |
| sa / lo / c4 BCF | 103.3 / 108.2 / 108.4% | 107.1 / 112.6 / 112.9% | 106.6 / 111.0 / 111.2% | 107.3 / 111.8 / 112.1% |
| cohort x.vcf.gz | 191.3% | 192.3% | 191.7% | 192.7% |

PGO recovers part of what the default layout loses (C against A: up to +3.3% on sa.bcf) but
does not beat the alignment flag, and on top of it is neutral to 2.6% worse (D against B).

**Ratio gate.** The 8:1 gate uses the output buffer's length; pasta passes exactly sized
buffers (the block's ISIZE), so the review's concern about oversized reusable buffers does
not apply to it.

**Outcome.** No change met the acceptance criterion; the crate is as before this round. By
the review's stopping rule, further decoder restructuring needs a newly measured
bottleneck first. On this CPU, the loop's cost per entry type is what every multi-symbol
idea since literal pairs has run into.

### 5.13 Apple silicon (M4 Pro): long reads, and a two-lane batch decoder

Everything above was measured on x86-64. On an M4 Pro (128 KiB L1D, 31 GPRs) the same
code was well short of its x86 lead (bgzipped VCF 110% against 170–191%), and nothing in
the codec was Arm-specific. Measured with paired, alternating runs of the profiling
example (median of 3; 24 MiB windows), on long-read BAMs: Revio HiFi whole genomes
(windows of two samples), an older HiFi TRGT BAM, and two ONT BAMs; four more windows
(two other HiFi samples, another window of the first, another ONT BAM) held out until
the end.

| Change | HiFi | ONT | |
|---|---|---|---|
| Baseline | 116% | 113–116% | |
| CRC32: PMULL folding plus the CRC unit at once (`crc32_arm_hybrid`) | — | — | kept: 90 vs 77 GiB/s at 64 KiB, but ~2% of decode time |
| 16-byte match copies (offset ≥ 16) | −1% | | reverted |
| `fastloop` out of line on Arm: inlined, LLVM kept the table pointer and output bound on the stack and reloaded them every symbol | +2–3 | +3 | kept |
| Literal pairs up to 16:1 instead of 8:1 on Arm (Illumina +3.5) | +0–1 | +0–1 | kept |
| Up to 4 literals per lookup (64-bit entries; simulated) | | | not built: 3–8% fewer lookups than pairs |
| Two blocks at once (`deflate_decompress_pair`) | 140% | 140–143% | superseded |
| The pair loop's bounds as a counted round budget (spills 10 → 8) | ±1 | ±1 | reverted |
| Three blocks at once (37 spills; budgeted bounds did not reduce them) | ±1 | +3–7 | reverted: TRGT 106% → 98% |
| **Two lanes over a batch** (`deflate_decompress_batch`), Arm only | **142–145%** | **145–150%** | **kept** |

**The batch decoder** answers §7 question 3. The fastloop's body runs as one step of
each of two streams per round, inlined twice so each stream's branches have their own
predictor sites; headers, tables and the generic loop run per stream (`Stream::advance`).
With fixed pairs, 18% (HiFi) to 24% (ONT) of the time was one stream running alone after
the shorter of the two ended; as a batch, a lane that ends takes the next stream, so the
lanes run alone only at the end of the batch. A stream that meets an error leaves its
lane and is decoded again alone, so errors are the single decoder's exactly. On x86-64 the
§5.6 version lost (87–99%) to spills and shared mispredicts; on the M4 the loop fits in
registers (657 instructions, 10 stack loads) and the reorder window overlaps the two
lookup chains. Batches of 8 to 256 blocks gave the same speed; 2 lost the gain on ONT.
With the CRC check, as pasta uses it, batches of 32:

| File | Single block | Batch |
|---|---|---|
| HiFi, tuning (2 windows) | 119.6% | 143.1–143.3% |
| HiFi, held out (3 windows) | 118.3–119.4% | 141.8–144.9% |
| ONT, tuning (2) | 114.7–118.8% | 145.0–150.3% |
| ONT, held out | 114.8% | 144.5% |
| HiFi TRGT (12:1) | 105.0% | 105.0% |

**Compression on long reads** (4 MiB windows, two HiFi samples and ONT): level 1 97–98%
(same output), level 5 132–146% (output 0.4–0.7% smaller on HiFi, equal on ONT), level 6
111–121%, level 9 118–120% (equal output). `prefetch` was a no-op on Arm, where
libdeflate's C issues `prfm`; issuing it (stable inline assembly, as `_prefetch` is not
stable) gave level 6 +1.5–4 points. Unchecked loads in level 1's matcher, as level 2–9's
has, gave under a point and were reverted.

**Validation**, on this code: the whole suite (29 tests: encoder round trips at every
level, damaged-stream agreement with libdeflate, guard pages, vectors), with batches
checked against single calls in batches of 1, 2, 3, 5, 7 and all, on damaged and
truncated streams, at buffer edges, and on every block of five long-read windows (8,080
blocks) — release, debug, portable and x86-64 (Rosetta) builds;
**AddressSanitizer** (nightly 2026-11-06, `-Zsanitizer=address -Zbuild-std`,
`detect_odr_violation=0` as in §3.5, 400 damage trials, 8 MiB windows of the five
files): no errors; **Miri** (Stacked Borrows) on `tests/vectors.rs`, its four tests now including
batches: no undefined behaviour (2 h 16 min; under Miri the CRC32 is slice-by-8 and the
encoder's prefetch is left out, as Miri runs neither PMULL/EOR3 nor inline assembly);
`cargo fmt` and `clippy -D warnings` clean on arm64 and
x86-64. On real x86 hardware, where the batch API decodes one by one: see the x86-64 check
below.

**Review round.** An in-depth review of this work found no decoding bugs; its fixes:

- The batch decoder's step was a copy of the fastloop's body. Both now run on one
  `Fast::step` (`fastloop` drives one stream, `fastloop_pair` two), and `fastloop_pair` is
  an `unsafe fn` whose contract is that both streams are ready. On the M4 the single path
  measured the same or slightly faster (HiFi 118–119% to 120–121%); batches 148–155% of
  libdeflate on long reads without the CRC. On x86-64 the BMI2 decoder grew from 897 to
  923 instructions, with a few more stack accesses in some loops: 2% slower there, mostly
  recovered since (the x86-64 check below).
- The hybrid CRC32's short inputs (under one 1,584-byte chunk) go to the PMULL code from
  the dispatch closure, a compare and a branch away: 19.8 GiB/s at 64 bytes against 18.8
  for the original code and 23.8 for libdeflate, which still leads at 64–128 bytes. The
  tuning generics became named constants.
- `deflate_decompress_batch` writes into a caller-owned `results` slice and keeps its
  scratch between calls, so a batch allocates nothing.
- The prefetch is declared `nomem` (no measurable change).
- A guard-page test no longer builds a one-byte-short output for empty blocks; build
  directories `target-*/` are ignored.

**x86-64 check** (Xeon Gold 5418Y, the decoder without the CRC; best of 5 interleaved runs
on an idle core, mean of three loop-alignment builds; % of libdeflate). The refactor
above cost x86-64 2%, on every file and every layout. The cause was register allocation,
not `Fast::step`: inlined into the decoder, the fastloop's registers were allocated with
the rest of the decoder's, and LLVM kept the table pointer on the stack and reloaded it
every symbol. (`main`'s loop body put back into this code was slower still: 0.833 against
0.848 GiB/s on Illumina.) With BMI2 the fastloop is now a function of its own,
`fastloop_bmi2`, as on Arm; Arm's code paths are unchanged.

| File | `main` | Refactor | `fastloop_bmi2` |
|---|---|---|---|
| GIAB HG001 ONT, window at 30% | 106.3% | 104.1% | 105.1% |
| HiFi chr20 | 107.2% | 105.6% | 106.6% |
| Arabidopsis HiFi (minimap2) | 125.1% | 122.7% | 125.2% |
| GIAB HG002 Revio, window at 25% | 107.9% | 105.8% | 106.9% |
| x.bam | 105.4% | 104.8% | 105.6% |
| Illumina chr20 | 107.4% | 104.9% | 106.8% |
| c3.bcf | 113.5% | 110.8% | 113.4% |
| x.bcf | 124.1% | 122.3% | 123.5% |
| 1kGP VCF | 169.9% | 169.0% | 170.0% |
| **Total (time-weighted)** | **113.1%** | **111.0%** | **112.4%** |

About a point is left on the GIAB long reads. The batch decoder's lanes, forced on for
the check, passed the batch tests in release and debug builds, and a libFuzzer target
(not in the repository) ran 2.47 million batches under AddressSanitizer (20 min, 4
workers), each stream checked against the single decoder, with fresh and reused
decompressors: no differences or errors.

### 5.14 Apple silicon check of `main` at `ff9b578` (level 1's new matcher)

`benches/arm-check.sh` on the M4 Pro (rustc 1.98.1), on `main` after the merge of §5.13
and the compressor commits that followed (fallible constructors, Rust 1.91, level 1's
short table and skipped inserts, level 0's stored-block fix), against the same check of
`c044c57` (§5.13's last commit). Three files: a PacBio HiFi BAM, an ONT BAM, and an
annotated VCF.

| % of libdeflate | HiFi BAM | ONT BAM | VCF |
|---|---|---|---|
| Level 1, `c044c57` | 96.0% | 96.8% | 97.7% |
| **Level 1, `ff9b578`** | **107.0%** | **121.1%** | **99.4%** |
| Level 1 size, `ff9b578` | 99.979% | 99.959% | 99.925% |
| Level 5 | 144.7% | 133.8% | 134.9% |
| Level 6 | 118.5% | 111.2% | 115.5% |
| Decompression with the CRC, one block per call | 123.2% | 119.8% | 112.9% |

Level 1 now beats libdeflate on long-read BAMs on Arm despite the 512 KiB short table
being four times L1; the VCF is even. Levels 5 and 6 and decompression are unchanged within
noise (§5.13's batch decoder is not in the script; it is unchanged). Level 5's VCF output
is 1.45% larger than libdeflate's, as before these commits. CRC32: 91.1 vs 78.0 GiB/s at
64 KiB, 71.3 vs 62.5 at 4 KiB, 21.4 vs 24.0 at 64 bytes.

Validation: the suite (34 tests) passes in release and debug builds, with Rust 1.91, as
x86-64 under Rosetta, and under AddressSanitizer (nightly, `-Zbuild-std`); 20,002 real
blocks decoded with no failures; `cargo fmt` and `clippy -D warnings` clean.

## 6. Where zflate's lead comes from, and its limits

- **Match-heavy data (x.bam-like BAM, BCF, VCF).** ~95% of decode time is Huffman decoding:
  each match is two dependent table lookups (literal/length, then offset), each roughly
  load → shift → mask → next load, ~7 cycles. libdeflate has the same chain. zflate's lead
  there comes from table building, which is a large share on BCF/VCF (hence 110–190%).
  On x.bam, table building and code-length reading are ~16% of the time: removing them
  *entirely* would give roughly 125%, halving them ~114%.
- **Literal-heavy data (real BAMs: long reads, Illumina).** The loop is most of the time,
  bounded by lookup latency and mispredicted literal/match branches (§5.11). Table building
  is 4% on ONT but ~15% on HiFi, whose BGZF blocks hold about three DEFLATE blocks each.
  Literal pairs help because they halve the lookups and branches on literal runs, not
  because they shorten the chain. The gain varies by data: most on Arabidopsis HiFi
  (125%) and the fly files (113%), least on the human GIAB ONT and HiFi files (106–108%).
  Pairs fit only where two literal codes total 12 bits or less, so shorter literal codes
  gain more.
- Every way tried of using the idle execution units — a second independent stream (twice:
  branchy, then branch-free), and the CRC — cost more than it gained on this CPU: the
  branchy versions in mispredicts and registers, the branch-free one in instructions
  (§5.11).
- **150% looks out of reach here**: on HiFi a token costs ~7 cycles of lookups, ~6–7 of
  mispredicts and ~3 of table building, against libdeflate's 20 in all; 150% would need
  13.4. Precode pairs, literal + length entries and fused entries built at their natural
  depth all lost (§5.11, §5.12): each new entry type slows the decode loop's common path by
  more than its saved lookups.
- Results are sensitive to code layout (±1.5% BAM, up to ±4% BCF) and stable Rust has no way
  to pin loop alignment.

## 7. Open questions

1. Literal pairs paid; literal + length entries and fused length + offset entries, tried
   as the second review proposed, did not (§5.12): each slowed the decode loop's common
   path by more than its saved lookups (3.7% and 17%), before builder costs. Is there a known
   way to add an entry kind without that cost on this kind of core, or is pairs-only the
   practical limit for a libdeflate-style loop?
2. ISA-L 2.31.1 is slower than zflate on every file, and slower than libdeflate on BAM and
   BCF except ONT (§5.12). Is there
   another decoder known to be faster than libdeflate on x86-64 that would show where the
   remaining headroom is?
3. Would the two-block interleaved decoder be expected to win on wider-register or
   different microarchitectures (ARM64, AMD Zen), and is it worth testing there?
   *Answered for ARM64 (§5.13): on an M4 Pro it wins, 142–150% of libdeflate on long
   reads; AMD Zen is untested.*
4. Are there layout-stabilisation techniques for Rust (function/loop alignment) usable on
   stable that would make results less build-dependent? Instrumentation PGO recovered part
   of the default layout's loss but did not beat `-align-all-functions=6` (§5.12).

## 8. Status and where things are

- The Apple silicon work (§5.13) is in `main`.
- Benchmarks: `benches/codec.rs`, `benches/crc32.rs`; profiling example
  `examples/profile.rs` (modes `zflate`, `libdeflate`, `zflate+crc`, `libdeflate+crc`,
  and `zflate-batch`, `zflate-batch+crc`, with `BATCH` blocks per call).
- The code of the reverted experiments, the measurement scripts and the ISA-L comparison
  harness are not published: they were written for one machine's paths.
- pasta integration: done. pasta decompresses, compresses and checks CRC32s with zflate,
  and its multithreaded BGZF reader hands blocks to the batch decoder on Arm.
