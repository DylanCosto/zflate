# Benchmarks

zflate against libdeflate 1.26, one core, on the BGZF blocks of public files: GIAB samples
(Illumina, PacBio HiFi and ONT BAMs), 1000 Genomes, cohort BCFs, and non-human genomes. Each
percentage is zflate's speed (or size) as a share of libdeflate's on the same blocks, in the
same process. [`PERFORMANCE-LOG.md`](PERFORMANCE-LOG.md) describes the files
and the methods (§2), and every change that was measured on the way here, including those
that did not work. A wider comparison, against other codecs on long-read BAM and CRAM, is in
[`LONG-READ-CODECS.md`](LONG-READ-CODECS.md). Commit hashes in these documents refer to zflate's development history before its public release; the public repository starts from that code as one commit.

The machines:

- **x86-64**: Intel Xeon Gold 5418Y (Sapphire Rapids, AVX-512, BMI2), a shared virtual
  machine; every run is pinned to an idle core.
- **Arm**: Apple M4 Pro.

Code layout alone moves a single build's results by 1–4%, so the main tables are the mean of
three builds with loops aligned to 16, 32 and 64 bytes ("layout-robust"); single-build
tables say so. The raw outputs of the 2026-11-07 compression runs are in
[`results/`](results/).

The code for `cargo bench` (criterion) is in [`../benches/`](../benches/); the comparisons
here come from `examples/profile.rs` (decompression), `examples/compress.rs` (compression),
`examples/validate.rs` (whole files) and `benches/arm-check.sh` (Apple silicon). All of them
take BGZF files as arguments, or from `ZFLATE_CORPUS`.

## Decompression

### x86-64, every corpus file

Layout-robust, 48 MB windows of whole blocks (from the performance log, §3.1). A later
change moved the BMI2 decoder loop out of line; nine of these files measured after it are
within 1.2 points of these figures (§5.13).

| File | Set | libdeflate GiB/s | zflate GiB/s | zflate / libdeflate |
|---|---|---:|---:|---:|
| GIAB HG001 ONT, window at 30% | long read, tuning | 0.372 | 0.396 | 106.5% |
| GIAB HG002 HiFi chr20 | long read, tuning | 0.757 | 0.811 | 107.2% |
| Arabidopsis HiFi (minimap2) | long read, tuning | 0.647 | 0.809 | 125.1% |
| GIAB HG001 ONT, window at 70% | long read, held out | 0.374 | 0.400 | 106.8% |
| GIAB HG002 HiFi whole genome, window at 25% | long read, held out | 0.890 | 0.960 | 107.9% |
| fly HiFi | long read, held out | 0.402 | 0.456 | 113.5% |
| fly dense reads | long read, held out | 0.403 | 0.455 | 112.8% |
| Arabidopsis HiFi, a second alignment | long read, held out | 0.648 | 0.814 | 125.6% |
| fly reads vs. targets | long read, held out | 0.370 | 0.419 | 113.2% |
| `x.bam` | BAM, tuning | 2.862 | 3.010 | 105.2% |
| GIAB Illumina chr20 | BAM, tuning | 0.802 | 0.864 | 107.7% |
| GIAB HG001 ONT, window at 50% | BAM, held out | 0.400 | 0.426 | 106.5% |
| Illumina, rewritten | BAM, held out | 0.786 | 0.843 | 107.3% |
| `short.name.bam` | BAM, held out | 3.204 | 3.480 | 108.6% |
| `c3.bcf` | BCF, tuning | 3.195 | 3.595 | 112.5% |
| cohort `x.bcf` | BCF, tuning | 2.585 | 3.193 | 123.5% |
| `hi.bcf` | BCF, tuning | 3.175 | 3.564 | 112.3% |
| cohort `la.bcf` | BCF, held out | 2.590 | 3.189 | 123.1% |
| cohort `sorted.bcf` | BCF, held out | 2.595 | 3.200 | 123.3% |
| `sa.bcf` | BCF, held out | 2.896 | 3.081 | 106.4% |
| `lo.bcf` | BCF, held out | 3.116 | 3.500 | 112.3% |
| `c4.bcf` | BCF, held out | 3.156 | 3.539 | 112.2% |
| 1000 Genomes chr22 VCF | VCF | 1.636 | 2.780 | 169.9% |
| cohort `x.vcf.gz` | VCF | 1.344 | 2.571 | 191.4% |

Over each set (time-weighted): long reads 111.3% tuning, 112.6% held out; BAM 107.1% and
112.0%; BCF 116.4% and 115.5%; VCF 181.1%. A bgzipped FASTA: 120.8%.

### x86-64, with the CRC32 check

Decompression followed by the block's CRC32, each library with its own (one build, best of 3):

| File | libdeflate GiB/s | zflate GiB/s | zflate / libdeflate |
|---|---:|---:|---:|
| GIAB HG001 ONT | 0.369 | 0.391 | 106.0% |
| GIAB HG002 HiFi chr20 | 0.747 | 0.791 | 105.9% |
| Arabidopsis HiFi | 0.640 | 0.788 | 123.1% |
| `x.bam` | 2.709 | 2.845 | 105.0% |
| GIAB Illumina chr20 | 0.791 | 0.843 | 106.6% |
| `c3.bcf` | 2.982 | 3.228 | 108.2% |
| cohort `x.bcf` | 2.459 | 2.954 | 120.1% |
| 1000 Genomes chr22 VCF | 1.571 | 2.601 | 165.6% |

### Whole files

Every block of each file decoded and checked against its CRC32 and length, then timed per
core with the CRC32 check, four workers at once (`examples/validate.rs`). Over the eleven
BAMs and two CRAMs checked (§3.3 lists them all), every one of the 16,420,083 BGZF blocks
and 682,890 CRAM gzip blocks was correct.

| File | Blocks | Compressed → inflated | zflate / libdeflate |
|---|---:|---|---:|
| GIAB Illumina HG003 NovaSeq 35×, whole genome | 4,852,927 | 46.80 → 315.83 GB | 104.4% |
| GIAB HG002 PacBio HiFi, whole genome | 5,351,415 | 67.27 → 267.52 GB | 113.1% |
| GIAB HG001 ONT, whole genome | 5,097,827 | 115.57 → 254.64 GB | 105.7% |
| GIAB HG002 PacBio HiFi chr20 | 120,041 | 1.52 → 6.00 GB | 112.6% |
| Arabidopsis HiFi, two alignments | 146,559 / 149,839 | 1.63 → 7.61 / 1.68 → 7.82 GB | 121.2% / 120.3% |
| CRAM 3.1, Illumina chr20: its gzip blocks | 19,246 | 0.06 → 0.30 GB | 119.2% |
| CRAM 3.0, ONT HG003, whole genome: its gzip blocks | 663,644 | 2.47 → 9.74 GB | 123.5% |

### Apple M4 Pro

Long-read BAMs, one block per call and through `deflate_decompress_batch`, which decodes two
blocks at a time on Arm (with the CRC32 check; §5.13):

| File | One block | Batch |
|---|---:|---:|
| HiFi, tuning (2 windows) | 119.6% | 143.1–143.3% |
| HiFi, held out (3 windows) | 118.3–119.4% | 141.8–144.9% |
| ONT, tuning (2 windows) | 114.7–118.8% | 145.0–150.3% |
| ONT, held out | 114.8% | 144.5% |
| HiFi, highly compressible (12:1) | 105.0% | 105.0% |

On a HiFi BAM, an ONT BAM and an annotated VCF, one block per call with the CRC32 check:
123.2%, 119.8% and 112.9% (§5.14).

## Compression

### Levels 1 and 5, x86-64

Layout-robust, 4 MiB windows; the size is zflate's output as a share of libdeflate's
([raw](results/x86-levels-1-and-5.txt)).

| File | Level 1 speed | Level 1 size | Level 5 speed | Level 5 size |
|---|---:|---:|---:|---:|
| GIAB Illumina chr20 BAM | 95.8% | 99.91% | 125.1% | 100.75% |
| `x.bam` | 92.2% | 99.97% | 155.6% | 99.54% |
| GIAB HG002 HiFi chr20 BAM | 100.6% | 99.97% | 124.5% | 101.12% |
| GIAB HG001 ONT BAM | 115.8% | 99.97% | 121.6% | 99.90% |
| Arabidopsis HiFi BAM | 117.6% | 100.00% | 126.7% | 100.11% |
| cohort `x.bcf` | 231.9% | 99.97% | 168.8% | 99.83% |
| `c3.bcf` | 224.2% | 99.97% | 165.1% | 101.43% |
| 1000 Genomes chr22 VCF | 142.1% | 99.97% | 171.9% | 99.11% |
| `idx.vcf.gz` | 143.4% | 99.96% | 176.5% | 101.41% |
| bgzipped FASTA | 88.7% | 99.97% | 134.8% | 99.01% |
| annotation TSV | 89.3% | 99.97% | 122.5% | 101.14% |

### Every level, x86-64

One build, 4 MiB windows ([raw](results/x86-every-level.txt)). Speed:

| Level | Illumina BAM | HiFi BAM | ONT BAM | BCF | VCF | FASTA |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 116% | 113% | 100% | 102% | 101% | 105% |
| 1 | 96% | 101% | 115% | 231% | 141% | 88% |
| 2 | 102% | 104% | 101% | 109% | 109% | 99% |
| 3 | 100% | 101% | 100% | 108% | 108% | 98% |
| 4 | 98% | 99% | 98% | 105% | 106% | 96% |
| 5 | 125% | 124% | 121% | 167% | 171% | 134% |
| 6 | 98% | 98% | 100% | 103% | 102% | 100% |
| 7 | 97% | 96% | 98% | 100% | 100% | 98% |
| 8 | 94% | 93% | 94% | 96% | 97% | 95% |
| 9 | 93% | 92% | 93% | 94% | 96% | 94% |
| 10 | 110% | 101% | 102% | 97% | 100% | 89% |
| 11 | 117% | 105% | 108% | 99% | 102% | 89% |
| 12 | 120% | 106% | 112% | 101% | 105% | 89% |

Size, the difference from libdeflate's output:

| Level | Illumina BAM | HiFi BAM | ONT BAM | BCF | VCF | FASTA |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 1 | −0.09% | −0.03% | −0.03% | −0.03% | −0.03% | −0.03% |
| 2–4 | 0.00% | 0.00% | 0.00% | −0.02 to −0.03% | −0.03 to −0.04% | 0.00% |
| 5 | +0.75% | +1.12% | −0.11% | −0.17% | −0.89% | −0.99% |
| 6–12 | 0.00% | 0.00% | 0.00 to +0.01% | −0.02 to −0.03% | −0.03 to −0.04% | 0.00% |

Level 5 is tuned for speed: it searches less deeply than libdeflate's and takes no 3-byte
matches. The other levels follow libdeflate's parameters.

### Apple M4 Pro

`benches/arm-check.sh`, 4 MiB windows (§5.14):

| | HiFi BAM | ONT BAM | VCF |
|---|---:|---:|---:|
| Level 1 | 107.0% | 121.1% | 99.4% |
| Level 1 size | 99.98% | 99.96% | 99.93% |
| Level 5 | 144.7% | 133.8% | 134.9% |
| Level 6 | 118.5% | 111.2% | 115.5% |

## CRC32

x86-64 (VPCLMULQDQ with AVX-512), GiB/s:

| Bytes | libdeflate | zflate | crc32fast |
|---:|---:|---:|---:|
| 16 | 3.7 | 3.8–4.0 | 2.7 |
| 64 | 12.2 | 10.4 | 7.6 |
| 256 | 34.1 | 33.1 | 20.1 |
| 1 KiB | 58.4 | 55.8 | 32.7 |
| 4 KiB | 71.2 | 70.2 | 55.1 |
| 64 KiB | 66.1 | 66.1 | 60.4 |
| 1 MiB | 76.6 | 66.7 | 60.7 |

Apple M4 Pro (PMULL with EOR3, and the CRC instructions at once), GiB/s: 91.1 against
libdeflate's 78.0 at 64 KiB, 71.3 against 62.5 at 4 KiB, and 21.4 against 24.0 at 64 bytes.

## Long-read codecs, x86-64

[Long-read codec comparison: level 5, x86](LONG-READ-CODECS.md) records the
2026-11-07 comparison across C and Rust codecs on ONT and HiFi BAM, CRAM 3.0 and
CRAM 3.1, including comparable-size settings and shared-input decoding. It includes
raw measurements, correctness exclusions and the limits of the comparison.
This is a separate single-build experiment from the layout-robust results above.
