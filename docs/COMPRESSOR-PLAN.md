# Compressor plan

zflate gets a DEFLATE compressor (raw, zlib, gzip) to replace libdeflate's in pasta. The output
need not match libdeflate's byte for byte. It must decompress to the input, with every
decoder, at every level, and be about as small as libdeflate's.

## What pasta needs

- BGZF blocks of at most 65,280 bytes, at libdeflate levels 0–12. HTSlib-style levels 1–9 map
  to libdeflate's {1, 2, 3, 5, 6, 7, 8, 10, 12}, so BAM and BCF are written at 7 by default,
  sort's temporary files at 1, and `CompressionLevel::default()` is 6.
- CRAM gzip blocks: gzip at levels 0–9, any size (often far larger than 64 KiB).
- A bound to size the output buffer. Data that does not compress must still fit in it, since
  the BGZF writer has no fallback with libdeflate: output is never larger than stored blocks.

## Design

libdeflate's, ported: its algorithms are well tuned, and matching its ratio is the bar.

| Levels | Match finder | Parser | Blocks |
|---|---|---|---|
| 0 | none | stored | 65,535-byte stored blocks |
| 1 | hash table, 2 entries a bucket | greedy, nice length 32 | at most 65,535 bytes or 8,192 sequences |
| 2–4 | hash chains (3- and 4-byte hashes) | greedy | split by statistics, at most 300,000 bytes |
| 5–7 | hash chains | lazy (one step) | as 2–4 |
| 8–9 | hash chains | lazy (two steps) | as 2–4 |
| 10–12 | binary trees | near-optimal, iterated costs | as 2–4 |

Shared by all levels:
- length-limited canonical Huffman codes (litlen 14 bits, offset 15, precode 7)
- the choice of the cheapest of dynamic, static and stored for each block, from its exact
  cost in bits; the cost also proves the block fits before it is written
- the block-split heuristic (literal and match statistics, checked every 512 observations)
- the minimum match length chosen from the literals' spread

Positions in the match finders are `i16` relative to a base that slides every 32 KiB, as in
libdeflate, so the tables stay small: 256 KiB for hash chains, 128 KiB for the hash table.
(Later, inputs up to 64 KiB kept absolute positions, and level 1's table for them grew to
512 KiB: see "Level 1, and inserts that cannot matter".)

Rust shape:
- `Compressor::new(CompressionLevel)`; `deflate_compress`, `zlib_compress`,
  `gzip_compress` into a caller's buffer, returning the length or
  `CompressionError::InsufficientSpace`; `*_compress_bound` functions
- the state per level is an enum of boxed match finders, allocated once and reused
- safe code first; `unsafe` only where a measurement shows it pays, each with its proof,
  as in the decoder

## Correctness

- unit tests for the Huffman builder (Kraft sum, length limit, canonical order), the precode
  and the slot tables
- round trips at every level through zflate, libdeflate and zlib: corpus blocks, and made-up
  inputs (empty, 1 byte, zeros, random, short and long periods, offsets at 32,768, inputs past
  32 KiB and 300 KiB, for the window slide and block splits)
- output never larger than the bound; too-small buffers give `InsufficientSpace`, never a
  short or broken stream
- Miri on small inputs, and the guard-page test, once there is `unsafe`

## Targets, per level, on the corpus

- size within 0.5% of libdeflate's
- speed at least libdeflate's, then better

## Stages

1. Levels 0–9, with tests and the comparison tool.
2. Levels 10–12.
3. Speed: profile against libdeflate.
4. pasta: compression through zflate with the `zflate` feature; tests; end-to-end timings.

## Status

Stages 1 and 2 are done: all 13 levels, raw DEFLATE, zlib and gzip. Output is 0.03–0.05%
smaller than libdeflate's at levels 1–12 on the test inputs (the first position is hashed
properly; libdeflate hashes it as 0), and identical in size at level 0.

Tests: round trips through zflate, libdeflate and zlib at every level and format; the bound;
`InsufficientSpace` one byte short, and output independent of the room given; the exact
block cost checked against the bits written (debug assertions); guard pages around inputs
and outputs; Miri on a zflate-only round trip.

### Speed experiments, levels 1–9

Each row is the zflate / libdeflate speed ratio on `ilmn.chr20.bam` and `x.bcf` (4 MB
windows), from the layout-robust method once it was adopted: the mean of three builds with
loops aligned to 16, 32 and 64 bytes (each build agreed within about 1%).

| Change | Effect |
|---|---|
| The safe port | 84–90% |
| Unchecked loads of chain candidates (clamped below `pos`) | +0–2 |
| Prefetching the next position's hash buckets, as libdeflate does | ±1 |
| Bit writer's slow path given values, so its state stays in registers | 0 |
| Sequences in a fixed slice, the counts in a local (no `Vec::push`) | +2 at level 1, 0 elsewhere |
| Each parser in its own function | 0 |
| Hash-chain tables in one allocation | −9 at levels 6–7: the link lookup took 4 instructions, not 2 |
| … with the links at offset 0 | back to 86–88% |
| A never-taken branch instead of the clamp | 0 |
| `-C target-cpu=native` (not kept) | 0 to −5 |
| `lz_extend` comparing 4 words unrolled, `cold_path` on libdeflate's unlikely branches | +3–5: 88–90% |
| Unchecked `lz_extend` for candidates (reverted) | −1 |
| Offset slot stored in the sequence, length tables padded, literal tail unrolled | ±2 |
| Block state as two pointers and a fixed-size array | ±1 |

The profile is flat (no source line above 7%), and both compressors spend the same share in
match finding and block writing, each about 10% slower in zflate. The loops spill to the
stack where libdeflate's do not, and changes that only move register allocation swing the
result by ±2–9%. The likely remaining gap is code generation: libdeflate here is built by
GCC 15, zflate by LLVM. Not tested: there is no clang on the machine to build libdeflate
with LLVM too.

### The VCF size bug

The minimum-match-length table (`choose_min_match_len`) was copied from libdeflate by hand
with two entries missing, so data with 43–44 distinct byte values got a minimum of 4
instead of 5: VCF output was 0.2–7% larger than libdeflate's at levels 4–9. The table is
now generated from libdeflate's source, and `sizes_track_libdeflate` checks each input
within 1%, not only the total, which had hidden it.

### Level 5

Level 5 is the target: libdeflate's level 5 keeps most of level 7's
ratio at about twice its speed. The goal: the same size, faster.

What does not work, measured against libdeflate level 5 on `ilmn.chr20.bam`,
`hifi.chr20.bam`, `x.bcf` and the 1000 Genomes chr22 VCF:

- **Other parameters.** A sweep of 79 combinations of parser (greedy, lazy, lazy2), depth
  and nice length found none as small as libdeflate's level 5 and faster: its settings
  are on the frontier. Even depth 2 is only ~1.2× faster, at +4–14% size, so most of level
  5's time is per position, not deep chain walks.
- **ISA-L** (igzip 2.31) is 1.6–5× faster but 12–30% larger than libdeflate level 5, larger
  even than libdeflate level 1.
- **A row match finder** (zstd's design: 32 recent positions and 8-bit tags a bucket):
  slower at levels 2–6 (inserting into 640 KB of tables) and larger at 6–9.
- **Paired chain walks** (two positions' chains interleaved): slower at level 5.
- **Inserting only part of each long match**: up to 1.6× faster on BCF, but +0.3–2% size.
  (A cap on every match; leaving out only the inserts no search can reach, in periodic
  matches, keeps the output identical: see "Level 1, and inserts that cannot matter".)
- Deferring the 3-byte candidate's load; huge pages for the tables; a BMI2 or native build;
  keeping the 3-byte heads only while the minimum match is 3 or 4 (VCF +7 points, but HiFi
  and BCF switch often and lose 4–8); smaller hash tables (VCF and BCF +2–3, BAM −2).

What works, all with identical output:

| Change | Effect (speed vs libdeflate level 5) |
|---|---|
| Start | 88–90% |
| Panic-free literal reads (`get(..).unwrap_or`) and a returning, not panicking, search guard | +2–9 points |
| Unchecked loads in the chain search and insertion, each proven by one check | +1–2 |
| Rebase with a saturating add (`paddsw`), as libdeflate's SSE2 version | +1 |
| Codeword reversal by byte table; a plain loop for precode runs | +1 |
| The block body written through a local `BitWriter` (kept in registers) | +0–0.5 |
| Inputs up to 64 KiB (BGZF's) keep absolute 16-bit positions: no window moves, no per-position window test, and one unsigned compare checks a candidate is in the window and before `pos` (replacing the cutoff test and the soundness guard) | +2–7 |
| Flush after a run's last literals only when there are some (runs between matches are often empty) | +0.2–0.6 |
| Codeword and length packed in one table entry for the writer | 0 (kept: simpler) |
| Block-split observations derived from the symbol counts at each check, not counted per symbol in memory (the same decisions) | +2–4 |
| Now | BAM 100.5–100.6%, BCF 105.6%, VCF 106.4% |

Also tried in this round, and reverted:

- Multi-hop chain links (libdeflate's PR 465, after 1.26: four next positions per link
  entry, a load every fourth step). As `[i16; 4]`, BCF fell to 45%: inserting a position
  in a run reads the entry the last insertion wrote, and four 2-byte stores do not forward
  to an 8-byte load. As `u64`, BCF 75%, BAM +1.5. With the run's link kept in a register,
  every file was slower (88–96%): the 256 KiB link table and the extra insertion work cost
  more than the shorter waits save here.
- Two tables of symbol counts, alternating, so that repeated symbols do not chain their
  increments (zstd's `HIST_count_parallel`): 1–1.5 points slower.

A survey of other compressors (zlib-ng, Chromium and Cloudflare zlib, AMD AOCL, Intel IPP,
zlib-rs, miniz_oxide, ISA-L, slz) found no published measurement of any beating libdeflate
at levels 5–6, per core, at equal size.

The lesson of the biggest one: a possible panic in a hot loop costs more than its compare;
it constrains the whole loop's code. Removing the remaining bounds checks altogether, as
an experiment, was worth 6–9 points, of which these changes recovered most.

The two guards left in the search (a candidate must be before `pos`; the search's bytes
must be in the input) cost 1.3–2.6 points. They are redundant with the tables' invariants,
but they make the unchecked loads sound whatever the tables hold. Without them: 95.6–99.2%.

### Level 5 trades a little size for speed

Level 5 may be up to 1–2% larger than libdeflate's for speed. libdeflate's level 5
searches 16 positions deep and its look-ahead 8, and stops at 30 bytes; zflate's searches 12
deep, its look-ahead 3, and stops at 48. Chosen from sweeps of depth, look-ahead depth and
nice length (look-ahead depth is the strongest lever, and libdeflate fixes it at half),
keeping every file within 1.5%.

### Level 5 takes no 3-byte matches

libdeflate allows 3-byte matches (within 8 KiB) when a block's first 4 KiB uses 80 or more
distinct byte values, which binary BAM and BCF data always do. There they cost more than
they save: with a minimum of 4, output at level 5 is 0.4% smaller on both ONT BAMs, 0.5%
on `x.bcf`, 0.1% on `c3.bcf`, and 0.15% and 0.08% larger on HiFi and Illumina BAMs. VCF,
TSV and FASTA already get a minimum of 4 or more. (A minimum of 5 is mixed: −1.9% on
`minimap2.bam`, but +1.8% on Illumina and +0.9% on HiFi.)

Without 3-byte matches, level 5 keeps no 3-byte hash table: no hash, head load and store at
every position, and a 64 KiB smaller table to clear for each input. That alone took the
speed from 107.6% to 120.9% on ONT, 113.2% to 121.1% on HiFi and 117% to 136% on BCF, at
the same search settings; VCF went from 119% to 137% with identical output.
Levels 2–4 and 6–9 keep libdeflate's minimum.

With libdeflate's own level 5 settings (16, 8, 30), this alone is 109–129% of libdeflate's
speed and within +0.14% of its size on every file (smaller on most); with (16, 8, 65),
smaller on every file at 106–123%.

Mean of three builds, against libdeflate 1.26 level 5:

| File | Speed | Size | Before this round |
|---|---|---|---|
| HiFi BAM (`hifi.chr20.bam`) | 123.7% | +1.12% | 114.0%, +0.79% |
| ONT BAM (HG001 PAW79146) | 123.0% | −0.11% | 108.0%, +0.24% |
| ONT BAM (`minimap2.bam`) | 126.4% | +0.11% | 112.5%, +0.40% |
| Illumina BAM (`ilmn.chr20.bam`) | 125.4% | +0.75% | 110.8%, +0.49% |
| BCF (`x.bcf`, `c3.bcf`) | 139.5%, 139.3% | −0.17%, +1.43% | 117.8%, 115.2% |
| VCF (1000 Genomes chr22, `idx.vcf.gz`) | 140.6%, 141.5% | −0.89%, +1.41% | 119.2%, 117.8% |
| Annotation TSV (`annot.tsv.gz`) | 122.2% | +1.14% | 115.6% |
| FASTA (`fragmented-reference.fa.gz`) | 133.9% | −0.99% | 129.4% |

The three builds agree within 1%, but for FASTA (129–138%). "Before this round": depth 12,
look-ahead 4, with 3-byte matches.

Tried in this round, and not kept:

- **Ending the look-ahead early.** Taking a match of 16, 24 or 32 bytes or more without
  looking ahead, or looking ahead less deeply after an 8-, 12- or 16-byte match (zlib's
  `max_lazy` and `good_match`), trades size for speed along the same curve as the depth
  settings: no better.
- **Counting symbols out of the search loop**, from the sequences, at each block-split check
  and at the block's end (identical output): 7% slower. Each literal run's length is
  unpredictable, so the separate loop mispredicts once a sequence, where the search loop
  had paid for that branch already; counting in four tables in turn did not help.
- **The lazy parser's threshold** (how much better the later match must be): 1 to 4
  makes no difference.
- **The block-split cutoff.** Splitting more readily (150, not 200) makes long-read BAMs
  0.25–0.45% smaller, but FASTA 1.0%, `c3.bcf` 0.3% and `idx.vcf.gz` 0.2% larger;
  splitting less (300) makes FASTA 2.5%, `c3.bcf` 1.2% and `idx.vcf.gz` 0.5% smaller, but
  long-read BAMs 0.65–0.9% larger. The best split depends on the data; a split test that
  estimates the bits saved, rather than comparing distributions, might win on both.

### Level 1, and inserts that cannot matter

Level 1 was the one level slower than libdeflate on most data (84–97% on BAM, VCF, FASTA
and TSV), and it writes pasta's temporary files. Two changes:

- **Inputs up to 64 KiB** (BGZF blocks) get their own level-1 table: 2^16 buckets of two
  `u32` entries, each an absolute 16-bit position below a 16-bit tag. The bucket index and
  the tag are the two halves of the 4-byte sequence times an odd constant, which is
  invertible, so a matching tag means the four bytes are equal: no load to check them, and
  the extension starts at byte 4. Longer inputs keep libdeflate's sliding `i16` table, in
  the first 128 KiB of the same allocation. Output is 0.00–0.09% smaller than before.
- **Inserts no search can reach** are left out, at levels 1 and 5 (the levels with no
  3-byte table). Inside a match of offset `p`, the sequence at each position is the one at
  `p` before it. Level 1 keeps two entries a bucket: if a match's last `2p` positions
  (before the three whose 4 bytes run past it) are inserted, every bucket an earlier
  position would have reached gets at least two newer entries, so the table ends as it
  would have. Level 5 searches 12 entries deep: its last `12p` positions leave every chain
  an earlier position joins with 12 newer entries in front of it, so no search reaches
  past them. Only matches of 31 or more positions take this path. A build with it turned
  off gives byte-identical output at every level, on the corpus and 400 synthetic inputs.

Mean of three builds (loops aligned to 16, 32, 64; best of three runs each; the builds
agree within 3.6 points), against libdeflate 1.26:

| File | Level 1 before | Level 1 after | Level 5 before | Level 5 after |
|---|---:|---:|---:|---:|
| `x.bcf`, `c3.bcf` | 158.6%, 159.9% | 231.9%, 224.2% | 138.8%, 138.2% | 168.8%, 165.1% |
| 1000 Genomes chr22 VCF, `idx.vcf.gz` | 89.3%, 89.9% | 142.1%, 143.4% | 140.2%, 141.0% | 171.9%, 176.5% |
| ONT BAM (HG001), `minimap2.bam` | 87.7%, 96.9% | 115.8%, 117.6% | 122.0%, 125.8% | 121.6%, 126.7% |
| HiFi BAM (`hifi.chr20.bam`) | 93.4% | 100.6% | 122.8% | 124.5% |
| Illumina BAM (`ilmn.chr20.bam`), `x.bam` | 88.4%, 83.7% | 95.8%, 92.2% | 125.1%, 155.9% | 125.1%, 155.6% |
| FASTA (`fragmented-reference.fa.gz`) | 80.2% | 88.7% | 134.7% | 134.8% |
| Annotation TSV (`annot.tsv.gz`) | 85.9% | 89.3% | 121.7% | 122.5% |

Level 5's output is unchanged; the gain is on data with long periodic runs, such as
genotype columns. Level 1 is still slower than libdeflate on short-read BAM, FASTA and TSV.

In pasta (the same commit built with each zflate): `vtools view -Ob` of a 500-sample BCF
1.13× faster on one thread and 1.07× with 8, the same bytes; a single-sample VCF.gz and
`bgzip` unchanged; `ntools sort` that spills to disk 1.02–1.04×; `markdup` 1.00–1.01×.

The costs: the short table is 512 KiB, four times libdeflate's. Clearing it for each input
is about 3% of level 1's time on `ilmn.chr20.bam` (a `memset`), and it no longer fits in
the 128 KiB L1 cache of Apple silicon's performance cores, as libdeflate's does. On an M4
Pro it is still faster: level 1 went from 96–98% of libdeflate to 107% (HiFi BAM), 121%
(ONT BAM) and 99.4% (VCF), with smaller output (`PERFORMANCE-LOG.md` §5.14).
