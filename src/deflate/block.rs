//! Ending a block: building its Huffman codes from the symbol frequencies, costing the three
//! block types exactly, and writing the cheapest.

use super::{
    bitstream::BitWriter,
    huffman::{codewords_from_lens, make_code},
    tables::{
        END_OF_BLOCK, FIRST_LEN_SYM, LENGTH_EXTRA_BITS, LENGTH_SLOT_BASE, LENGTH_TABLE_LEN,
        MAX_LITLEN_CODEWORD_LEN, MAX_MATCH_LEN, MAX_OFFSET_CODEWORD_LEN, MAX_PRE_CODEWORD_LEN,
        MIN_MATCH_LEN, NUM_LITERALS, NUM_LITLEN_SYMS, NUM_OFFSET_SLOTS, NUM_OFFSET_SYMS,
        NUM_PRECODE_SYMS, OFFSET_EXTRA_BITS, OFFSET_SLOT_BASE, OFFSET_SLOT_TABLE_LEN,
        PRECODE_EXTRA_BITS, PRECODE_LENS_PERMUTATION, length_slot, offset_slot, static_litlen_len,
    },
};

/// A run of literals and the match after it. The block's last sequence has no match
/// (length 0): it holds the literals after the last match.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Sequence {
    pub(super) litrunlen: u32,
    /// The match's length, and its offset's slot above bit 9, so that the writer need not
    /// find the slot again; 8 bytes a sequence, as libdeflate's.
    length_and_slot: u16,
    offset: u16,
}

impl Sequence {
    pub(super) fn new(litrunlen: u32, length: u32, offset: u32, offset_slot: usize) -> Self {
        debug_assert!(length <= MAX_MATCH_LEN && offset <= 32768 && offset_slot < 32);
        Self {
            litrunlen,
            length_and_slot: (length | (offset_slot as u32) << 9) as u16,
            offset: offset as u16,
        }
    }

    /// The literals after a block's last match.
    pub(super) fn last(litrunlen: u32) -> Self {
        Self::new(litrunlen, 0, 0, 0)
    }

    fn length(self) -> u32 {
        u32::from(self.length_and_slot) & 0x1ff
    }

    fn offset_slot(self) -> usize {
        usize::from(self.length_and_slot >> 9)
    }
}

/// How often each symbol occurs in the block so far.
#[derive(Clone)]
pub(super) struct Freqs {
    pub(super) litlen: [u32; NUM_LITLEN_SYMS],
    pub(super) offset: [u32; NUM_OFFSET_SYMS],
}

impl Freqs {
    pub(super) fn reset(&mut self) {
        self.litlen.fill(0);
        self.offset.fill(0);
    }

    #[inline(always)]
    pub(super) fn add_literal(&mut self, lit: u8) {
        self.litlen[usize::from(lit)] += 1;
    }

    #[inline(always)]
    pub(super) fn add_match(&mut self, length: u32, offset_slot: usize) {
        self.litlen[FIRST_LEN_SYM + length_slot(length)] += 1;
        self.offset[offset_slot & 31] += 1;
    }
}

/// A block's literal/length and offset codes.
#[derive(Clone)]
struct Codes {
    litlen_lens: [u8; NUM_LITLEN_SYMS],
    offset_lens: [u8; NUM_OFFSET_SYMS],
    litlen_codewords: [u32; NUM_LITLEN_SYMS],
    offset_codewords: [u32; NUM_OFFSET_SYMS],
}

impl Codes {
    const fn zeroed() -> Self {
        Self {
            litlen_lens: [0; NUM_LITLEN_SYMS],
            offset_lens: [0; NUM_OFFSET_SYMS],
            litlen_codewords: [0; NUM_LITLEN_SYMS],
            offset_codewords: [0; NUM_OFFSET_SYMS],
        }
    }

    /// RFC 1951's fixed codes.
    fn fixed() -> Self {
        let mut codes = Self::zeroed();
        for (sym, len) in codes.litlen_lens.iter_mut().enumerate() {
            *len = static_litlen_len(sym);
        }
        codes.offset_lens.fill(5);
        codewords_from_lens(&codes.litlen_lens, &mut codes.litlen_codewords);
        codewords_from_lens(&codes.offset_lens, &mut codes.offset_codewords);
        codes
    }
}

/// A dynamic block's code lengths, run-length coded with the precode.
struct Precode {
    freqs: [u32; NUM_PRECODE_SYMS],
    lens: [u8; NUM_PRECODE_SYMS],
    codewords: [u32; NUM_PRECODE_SYMS],
    /// Each a precode symbol, with its extra bits' value above bit 5.
    items: [u32; NUM_LITLEN_SYMS + NUM_OFFSET_SYMS],
    num_items: usize,
    num_litlen_syms: usize,
    num_offset_syms: usize,
    /// How many of the precode's lengths the header gives, in permuted order.
    num_explicit_lens: usize,
}

/// A position's choice on the near-optimal parser's path: a literal (`len` 1, the literal
/// above bit 9) or a match (`len`, the offset above bit 9), and the cost from there to the
/// block's end.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct OptimumNode {
    pub(super) cost_to_end: u32,
    pub(super) item: u32,
}

pub(super) const OPTIMUM_OFFSET_SHIFT: u32 = 9;
pub(super) const OPTIMUM_LEN_MASK: u32 = (1 << OPTIMUM_OFFSET_SHIFT) - 1;

/// What a block is made of: sequences of literals and matches, from the greedy and lazy
/// parsers, or the near-optimal parser's path through the block's nodes.
pub(super) enum Body<'a> {
    Sequences(&'a [Sequence]),
    Optimum(&'a [OptimumNode]),
}

/// The state for ending blocks, kept between them and between calls.
pub(super) struct BlockCoder {
    pub(super) freqs: Freqs,
    codes: Codes,
    fixed: Codes,
    precode: Precode,
    /// The block's codewords as the writer wants them.
    write: WriteTables,
}

/// Each symbol's codeword with its length above bit 24 ([`packed`]): one load a symbol
/// for the writer, not two. A match length's entry has its extra bits after its codeword.
struct WriteTables {
    literal: [u32; NUM_LITERALS],
    length: [u32; LENGTH_TABLE_LEN],
    offset: [u32; OFFSET_SLOT_TABLE_LEN],
    end_of_block: u32,
}

/// `bits` (at most 24 of them) and their count, in one word.
const fn packed(bits: u32, len: u32) -> u32 {
    bits | len << 24
}

/// The bits and their count from a [`packed`] entry.
#[inline(always)]
fn unpacked(entry: u32) -> (u32, u32) {
    (entry & 0xff_ffff, entry >> 24)
}

impl WriteTables {
    const fn zeroed() -> Self {
        Self {
            literal: [0; NUM_LITERALS],
            length: [0; LENGTH_TABLE_LEN],
            offset: [0; OFFSET_SLOT_TABLE_LEN],
            end_of_block: 0,
        }
    }

    fn fill(&mut self, codes: &Codes) {
        let litlen = |sym: usize| {
            (
                codes.litlen_codewords[sym],
                u32::from(codes.litlen_lens[sym]),
            )
        };
        for (lit, entry) in self.literal.iter_mut().enumerate() {
            let (codeword, len) = litlen(lit);
            *entry = packed(codeword, len);
        }
        for len in MIN_MATCH_LEN..=MAX_MATCH_LEN {
            let slot = length_slot(len);
            let (codeword, codeword_len) = litlen(FIRST_LEN_SYM + slot);
            let extra = len - LENGTH_SLOT_BASE[slot];
            self.length[len as usize] = packed(
                codeword | extra << codeword_len,
                codeword_len + u32::from(LENGTH_EXTRA_BITS[slot]),
            );
        }
        for (slot, entry) in self.offset.iter_mut().enumerate().take(NUM_OFFSET_SYMS) {
            *entry = packed(
                codes.offset_codewords[slot],
                u32::from(codes.offset_lens[slot]),
            );
        }
        let (codeword, len) = litlen(END_OF_BLOCK);
        self.end_of_block = packed(codeword, len);
    }
}

impl BlockCoder {
    pub(super) fn new() -> Self {
        Self {
            freqs: Freqs {
                litlen: [0; NUM_LITLEN_SYMS],
                offset: [0; NUM_OFFSET_SYMS],
            },
            codes: Codes::zeroed(),
            fixed: Codes::fixed(),
            precode: Precode {
                freqs: [0; NUM_PRECODE_SYMS],
                lens: [0; NUM_PRECODE_SYMS],
                codewords: [0; NUM_PRECODE_SYMS],
                items: [0; NUM_LITLEN_SYMS + NUM_OFFSET_SYMS],
                num_items: 0,
                num_litlen_syms: 0,
                num_offset_syms: 0,
                num_explicit_lens: 0,
            },
            write: WriteTables::zeroed(),
        }
    }

    /// Ends the block `data`, whose literals and matches `body` gives and whose symbols
    /// `self.freqs` counts: adds the end-of-block symbol, builds the codes, and writes the
    /// block as whichever type is smallest.
    pub(super) fn finish_block(
        &mut self,
        w: &mut BitWriter<'_>,
        data: &[u8],
        body: Body<'_>,
        is_final: bool,
    ) {
        self.freqs.litlen[END_OF_BLOCK] += 1;
        self.make_codes();
        self.flush_block(w, data, body, is_final);
    }

    /// Builds the block's codes from `self.freqs`, which counts the end-of-block symbol.
    pub(super) fn make_codes(&mut self) {
        let codes = &mut self.codes;
        make_code(
            MAX_LITLEN_CODEWORD_LEN,
            &self.freqs.litlen,
            &mut codes.litlen_lens,
            &mut codes.litlen_codewords,
        );
        make_code(
            MAX_OFFSET_CODEWORD_LEN,
            &self.freqs.offset,
            &mut codes.offset_lens,
            &mut codes.offset_codewords,
        );
    }

    /// Counts `block` as literals only, and builds codes for that.
    pub(super) fn choose_all_literals(&mut self, block: &[u8]) {
        self.freqs.reset();
        for &lit in block {
            self.freqs.add_literal(lit);
        }
        self.freqs.litlen[END_OF_BLOCK] += 1;
        self.make_codes();
    }

    /// The dynamic block's cost in bits with the codes built, but for BFINAL and BTYPE: what
    /// the near-optimal parser minimizes.
    pub(super) fn true_cost(&mut self) -> u32 {
        self.build_precode();
        self.costs(0, 0).0 - 3
    }

    /// The code lengths built, literal/length and offset.
    pub(super) fn lens(&self) -> (&[u8; NUM_LITLEN_SYMS], &[u8; NUM_OFFSET_SYMS]) {
        (&self.codes.litlen_lens, &self.codes.offset_lens)
    }

    /// The static codes' lengths.
    pub(super) fn fixed_lens(&self) -> (&[u8; NUM_LITLEN_SYMS], &[u8; NUM_OFFSET_SYMS]) {
        (&self.fixed.litlen_lens, &self.fixed.offset_lens)
    }

    /// Writes the block with the codes built, choosing its type by cost.
    pub(super) fn flush_block(
        &mut self,
        w: &mut BitWriter<'_>,
        data: &[u8],
        body: Body<'_>,
        is_final: bool,
    ) {
        debug_assert!(w.bitcount() <= 7 && !w.overflowed());
        self.build_precode();
        let (dynamic, fixed, stored) = self.costs(data.len(), w.bitcount());
        let best = dynamic.min(fixed).min(stored);
        // Knowing the block's size, check once that it fits.
        if (w.bitcount() as usize + best as usize).div_ceil(8) > w.room() {
            w.set_overflow();
            return;
        }
        let start = w.bits_written();
        if best == stored {
            // DEFLATE's stored blocks hold at most 65,535 bytes.
            let mut pieces = data.chunks(usize::from(u16::MAX)).peekable();
            while let Some(piece) = pieces.next() {
                w.write_stored(is_final && pieces.peek().is_none(), piece);
            }
        } else {
            let use_fixed = best == fixed;
            if use_fixed {
                w.add(u32::from(is_final), 1);
                w.add(1, 2);
                w.flush();
            } else {
                self.write_dynamic_header(w, is_final);
            }
            let codes = if use_fixed { &self.fixed } else { &self.codes };
            self.write.fill(codes);
            // The body's loop works on a local writer, which the compiler keeps in
            // registers; through `w`, it would store the buffer back after every sequence.
            let mut local = std::mem::replace(w, BitWriter::new(&mut []));
            match body {
                Body::Sequences(seqs) => write_sequences(&mut local, &self.write, data, seqs),
                Body::Optimum(nodes) => write_optimum(&mut local, &self.write, data.len(), nodes),
            }
            *w = local;
            let (bits, len) = unpacked(self.write.end_of_block);
            w.add(bits, len);
            w.flush();
        }
        debug_assert!(
            w.overflowed() || w.bits_written() - start == best as usize,
            "the block took {} bits, not {best}",
            w.bits_written() - start
        );
    }

    /// The cost in bits of the block as dynamic, fixed (static) and stored blocks, from the
    /// bit position `bitcount` within a byte.
    fn costs(&self, block_len: usize, bitcount: u32) -> (u32, u32, u32) {
        let (freqs, codes, fixed, pre) = (&self.freqs, &self.codes, &self.fixed, &self.precode);
        // Every block starts with BFINAL and BTYPE.
        let mut dynamic = 3 + 5 + 5 + 4 + 3 * pre.num_explicit_lens as u32;
        let mut fixed_cost = 3;
        for ((&freq, &len), &extra) in pre.freqs.iter().zip(&pre.lens).zip(&PRECODE_EXTRA_BITS) {
            dynamic += freq * (u32::from(extra) + u32::from(len));
        }
        // Literals and the end-of-block symbol, which has no extra bits.
        for sym in 0..=END_OF_BLOCK {
            let freq = if sym == END_OF_BLOCK {
                1
            } else {
                freqs.litlen[sym]
            };
            dynamic += freq * u32::from(codes.litlen_lens[sym]);
            fixed_cost += freq * u32::from(fixed.litlen_lens[sym]);
        }
        for (slot, &extra) in LENGTH_EXTRA_BITS.iter().enumerate() {
            let sym = FIRST_LEN_SYM + slot;
            let extra = u32::from(extra);
            dynamic += freqs.litlen[sym] * (extra + u32::from(codes.litlen_lens[sym]));
            fixed_cost += freqs.litlen[sym] * (extra + u32::from(fixed.litlen_lens[sym]));
        }
        for (slot, &extra) in OFFSET_EXTRA_BITS[..NUM_OFFSET_SLOTS].iter().enumerate() {
            let extra = u32::from(extra);
            dynamic += freqs.offset[slot] * (extra + u32::from(codes.offset_lens[slot]));
            fixed_cost += freqs.offset[slot] * (extra + u32::from(fixed.offset_lens[slot]));
        }
        // A stored block pads to a byte after its header, then has LEN and NLEN; each
        // further 65,535 bytes take another 5-byte header.
        let pieces = block_len.div_ceil(usize::from(u16::MAX)).max(1) as u32;
        let stored =
            3 + ((bitcount + 3).wrapping_neg() & 7) + 32 + 40 * (pieces - 1) + 8 * block_len as u32;
        (dynamic, fixed_cost, stored)
    }

    /// Run-length codes the dynamic code's lengths and builds the precode for them.
    fn build_precode(&mut self) {
        let codes = &self.codes;
        let pre = &mut self.precode;
        // Trailing unused symbols need not be sent, down to 257 literal/length symbols and 1
        // offset symbol.
        pre.num_litlen_syms = 257
            + codes.litlen_lens[257..]
                .iter()
                .rposition(|&l| l != 0)
                .map_or(0, |i| i + 1);
        pre.num_offset_syms = 1 + codes.offset_lens[1..]
            .iter()
            .rposition(|&l| l != 0)
            .map_or(0, |i| i + 1);
        // The two codes' lengths are run-length coded as one sequence.
        let mut lens = [0u8; NUM_LITLEN_SYMS + NUM_OFFSET_SYMS];
        let num_lens = pre.num_litlen_syms + pre.num_offset_syms;
        lens[..pre.num_litlen_syms].copy_from_slice(&codes.litlen_lens[..pre.num_litlen_syms]);
        lens[pre.num_litlen_syms..num_lens]
            .copy_from_slice(&codes.offset_lens[..pre.num_offset_syms]);
        pre.num_items = precode_items(&lens[..num_lens], &mut pre.freqs, &mut pre.items);
        make_code(
            MAX_PRE_CODEWORD_LEN,
            &pre.freqs,
            &mut pre.lens,
            &mut pre.codewords,
        );
        // Trailing zero lengths, in the header's order, need not be sent either, down to 4.
        pre.num_explicit_lens = 4 + PRECODE_LENS_PERMUTATION[4..]
            .iter()
            .rposition(|&sym| pre.lens[usize::from(sym)] != 0)
            .map_or(0, |i| i + 1);
    }

    /// Writes a dynamic block's header: the code sizes, the precode's lengths, and the
    /// code lengths in the precode.
    fn write_dynamic_header(&self, w: &mut BitWriter<'_>, is_final: bool) {
        let pre = &self.precode;
        w.add(u32::from(is_final), 1);
        w.add(2, 2);
        w.add(pre.num_litlen_syms as u32 - 257, 5);
        w.add(pre.num_offset_syms as u32 - 1, 5);
        w.add(pre.num_explicit_lens as u32 - 4, 4);
        // With the 7 bits that may be pending, 18 lengths more fit in the buffer: write the
        // first now, and the rest after a flush.
        let precode_len = |i: usize| u32::from(pre.lens[usize::from(PRECODE_LENS_PERMUTATION[i])]);
        w.add(precode_len(0), 3);
        w.flush();
        for i in 1..pre.num_explicit_lens {
            w.add(precode_len(i), 3);
        }
        w.flush();
        for &item in &pre.items[..pre.num_items] {
            let sym = (item & 0x1f) as usize;
            w.add(pre.codewords[sym], u32::from(pre.lens[sym]));
            w.add(item >> 5, u32::from(PRECODE_EXTRA_BITS[sym]));
            w.flush();
        }
    }
}

/// Run-length codes `lens` with the precode's symbols: 16 repeats the previous length 3 to
/// 6 times, 17 gives 3 to 10 zeros and 18 gives 11 to 138. Fills `freqs` and `items`, and
/// returns the number of items.
fn precode_items(lens: &[u8], freqs: &mut [u32; NUM_PRECODE_SYMS], items: &mut [u32]) -> usize {
    freqs.fill(0);
    let mut n = 0;
    let mut push = |item: u32, freqs: &mut [u32; NUM_PRECODE_SYMS]| {
        freqs[(item & 0x1f) as usize] += 1;
        items[n] = item;
        n += 1;
    };
    let mut start = 0;
    while start < lens.len() {
        let len = lens[start];
        let mut end = start + 1;
        while end < lens.len() && lens[end] == len {
            end += 1;
        }
        if len == 0 {
            while end - start >= 11 {
                let extra = (end - start - 11).min(0x7f);
                push(18 | (extra as u32) << 5, freqs);
                start += 11 + extra;
            }
            if end - start >= 3 {
                let extra = (end - start - 3).min(0x7);
                push(17 | (extra as u32) << 5, freqs);
                start += 3 + extra;
            }
        } else if end - start >= 4 {
            push(u32::from(len), freqs);
            start += 1;
            while end - start >= 3 {
                let extra = (end - start - 3).min(0x3);
                push(16 | (extra as u32) << 5, freqs);
                start += 3 + extra;
            }
        }
        for _ in start..end {
            push(u32::from(len), freqs);
        }
        start = end;
    }
    n
}

/// Writes the block's literals and matches.
fn write_sequences(w: &mut BitWriter<'_>, tables: &WriteTables, data: &[u8], seqs: &[Sequence]) {
    let literal = |w: &mut BitWriter<'_>, lit: u8| {
        let (bits, len) = unpacked(tables.literal[usize::from(lit)]);
        w.add(bits, len);
    };
    let mut pos = 0;
    for seq in seqs {
        let run = &data[pos..pos + seq.litrunlen as usize];
        pos += run.len();
        // Four literals of at most 14 bits fit with the 7 pending bits.
        let (quads, rest) = run.as_chunks::<4>();
        for quad in quads {
            for &lit in quad {
                literal(w, lit);
            }
            w.flush();
        }
        // The last one to three, without a loop: runs are mostly short, and often empty,
        // between matches, which needs no flush.
        if let [a, ref more @ ..] = *rest {
            literal(w, a);
            if let [b, ref more @ ..] = *more {
                literal(w, b);
                if let [c, ..] = *more {
                    literal(w, c);
                }
            }
            w.flush();
        }
        let length = seq.length();
        if length == 0 {
            break;
        }
        write_match(w, tables, length, u32::from(seq.offset), seq.offset_slot());
        pos += length as usize;
    }
    debug_assert_eq!(pos, data.len());
}

/// Writes the near-optimal parser's path through the block's `len` bytes.
fn write_optimum(w: &mut BitWriter<'_>, tables: &WriteTables, len: usize, nodes: &[OptimumNode]) {
    let mut pos = 0;
    while pos < len {
        let item = nodes[pos].item;
        let (length, offset) = (item & OPTIMUM_LEN_MASK, item >> OPTIMUM_OFFSET_SHIFT);
        if length == 1 {
            let (bits, n) = unpacked(tables.literal[offset as usize & 0xff]);
            w.add(bits, n);
            w.flush();
        } else {
            write_match(w, tables, length, offset, offset_slot(offset));
        }
        pos += length as usize;
    }
    debug_assert_eq!(pos, len);
}

/// Writes a match: its length's codeword and extra bits (at most 19 bits), its offset's
/// codeword (15) and extra bits (13), 47 bits with the 7 that may be pending.
#[inline(always)]
fn write_match(w: &mut BitWriter<'_>, tables: &WriteTables, length: u32, offset: u32, slot: usize) {
    let (bits, len) = unpacked(tables.length[length as usize & (LENGTH_TABLE_LEN - 1)]);
    w.add(bits, len);
    let slot = slot & 31;
    let (bits, len) = unpacked(tables.offset[slot]);
    w.add(bits, len);
    w.add(
        offset - OFFSET_SLOT_BASE[slot],
        u32::from(OFFSET_EXTRA_BITS[slot]),
    );
    w.flush();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn precode_runs() {
        let mut freqs = [0; NUM_PRECODE_SYMS];
        let mut items = [0; NUM_LITLEN_SYMS + NUM_OFFSET_SYMS];
        // 5 eights: one literal, then 16 repeating it 4 times.
        let n = precode_items(&[8; 5], &mut freqs, &mut items);
        assert_eq!(items[..n], [8, 16 | 1 << 5]);
        // 3 eights stay literal; 2 zeros stay literal; 140 zeros: 18 (138), then 2 literal.
        let mut lens = vec![8, 8, 8, 0, 0, 7];
        lens.extend([0; 140]);
        let n = precode_items(&lens, &mut freqs, &mut items);
        assert_eq!(items[..n], [8, 8, 8, 0, 0, 7, 18 | 127 << 5, 0, 0]);
        // 13 zeros: one 18.
        let n = precode_items(&[0; 13], &mut freqs, &mut items);
        assert_eq!(items[..n], [18 | 2 << 5]);
        // 10 zeros: 17 (10).
        let n = precode_items(&[0; 10], &mut freqs, &mut items);
        assert_eq!(items[..n], [17 | 7 << 5]);
        assert_eq!(freqs[17], 1);
    }
}
