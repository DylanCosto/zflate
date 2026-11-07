//! Huffman decode tables, in the format of libdeflate's `deflate_decompress.c` (MIT).
//!
//! A table maps the next `table_bits` bits of input, low bit first, to an entry for the
//! symbol whose codeword is a prefix of them. Codewords longer than `table_bits` go through
//! one level of subtables. Each entry carries everything decoding needs, so a symbol costs
//! one lookup:
//!
//! - bits 0-7: the bits to consume (the codeword length, plus the extra bits for lengths
//!   and offsets); for a subtable pointer, the main table bits.
//! - bits 8-11: the codeword length alone (to find the extra bits after it); for a subtable
//!   pointer, the subtable's bits.
//! - flags: [`LITERAL`] in bit 31; [`EXCEPTIONAL`] for the rare cases, with
//!   [`SUBTABLE_POINTER`] or [`END_OF_BLOCK`].
//! - bits 16 and up: the literal, the length or offset base, the precode symbol, or the
//!   subtable's start.

/// A literal: bits 16-23 are its value.
pub(super) const LITERAL: u32 = 0x8000_0000;
/// A subtable pointer or the end of the block.
pub(super) const EXCEPTIONAL: u32 = 0x8000;
pub(super) const SUBTABLE_POINTER: u32 = 0x4000;
pub(super) const END_OF_BLOCK: u32 = 0x2000;
/// With [`LITERAL`]: a second literal, in bits 8-15, follows the first in the entry's bits;
/// bits 0-7 are both codes' bits. Bits 8-15 then hold no flags.
pub(super) const PAIR: u32 = 0x4000_0000;

pub(super) const PRECODE_TABLEBITS: u32 = 7;
/// The most bits a literal/length main table has: 12 where it holds literal pairs, which
/// then fit more often; [`LITLEN_TABLEBITS_PLAIN`] otherwise.
pub(super) const LITLEN_TABLEBITS: u32 = 12;
/// The literal/length main table's bits without pairs, as libdeflate's.
pub(super) const LITLEN_TABLEBITS_PLAIN: u32 = 11;
pub(super) const OFFSET_TABLEBITS: u32 = 8;

/// The most entries each table can need, subtables included (zlib's `enough` program, as
/// libdeflate computes them).
pub(super) const PRECODE_ENOUGH: usize = 128;
// For 12 bits, a bound rather than zlib's exact `enough`: the main table, and at most one
// subtable of at most 2^(15 - 12) entries per symbol.
pub(super) const LITLEN_ENOUGH: usize =
    (1 << LITLEN_TABLEBITS) + NUM_LITLEN_SYMS * (1 << (15 - LITLEN_TABLEBITS));
pub(super) const OFFSET_ENOUGH: usize = 402;

pub(super) const NUM_PRECODE_SYMS: usize = 19;
pub(super) const NUM_LITLEN_SYMS: usize = 288;
pub(super) const NUM_OFFSET_SYMS: usize = 32;

pub(super) const MAX_PRE_CODEWORD_LEN: u32 = 7;
pub(super) const MAX_CODEWORD_LEN: u32 = 15;

/// The static part of each precode symbol's entry.
pub(super) static PRECODE_RESULTS: [u32; NUM_PRECODE_SYMS] = {
    let mut r = [0; NUM_PRECODE_SYMS];
    let mut sym = 0;
    while sym < NUM_PRECODE_SYMS {
        r[sym] = (sym as u32) << 16;
        sym += 1;
    }
    r
};

/// The static part of each literal/length symbol's entry: literals, end of block, then the
/// lengths' bases and extra bits (RFC 1951 3.2.5). Symbols 286 and 287 are not in the RFC;
/// like libdeflate, they decode as length 258.
pub(super) static LITLEN_RESULTS: [u32; NUM_LITLEN_SYMS] = {
    #[rustfmt::skip]
    const LENGTHS: [(u32, u32); 31] = [
        (3, 0), (4, 0), (5, 0), (6, 0), (7, 0), (8, 0), (9, 0), (10, 0),
        (11, 1), (13, 1), (15, 1), (17, 1), (19, 2), (23, 2), (27, 2), (31, 2),
        (35, 3), (43, 3), (51, 3), (59, 3), (67, 4), (83, 4), (99, 4), (115, 4),
        (131, 5), (163, 5), (195, 5), (227, 5), (258, 0), (258, 0), (258, 0),
    ];
    let mut r = [0; NUM_LITLEN_SYMS];
    let mut sym = 0;
    while sym < 256 {
        r[sym] = LITERAL | (sym as u32) << 16;
        sym += 1;
    }
    r[256] = EXCEPTIONAL | END_OF_BLOCK;
    let mut i = 0;
    while i < LENGTHS.len() {
        r[257 + i] = LENGTHS[i].0 << 16 | LENGTHS[i].1;
        i += 1;
    }
    r
};

/// The static part of each offset symbol's entry: bases and extra bits (RFC 1951 3.2.5).
/// Symbols 30 and 31 are not in the RFC; like libdeflate, they decode as symbol 29.
pub(super) static OFFSET_RESULTS: [u32; NUM_OFFSET_SYMS] = {
    #[rustfmt::skip]
    const OFFSETS: [(u32, u32); NUM_OFFSET_SYMS] = [
        (1, 0), (2, 0), (3, 0), (4, 0), (5, 1), (7, 1), (9, 2), (13, 2),
        (17, 3), (25, 3), (33, 4), (49, 4), (65, 5), (97, 5), (129, 6), (193, 6),
        (257, 7), (385, 7), (513, 8), (769, 8), (1025, 9), (1537, 9), (2049, 10), (3073, 10),
        (4097, 11), (6145, 11), (8193, 12), (12289, 12), (16385, 13), (24577, 13),
        (24577, 13), (24577, 13),
    ];
    let mut r = [0; NUM_OFFSET_SYMS];
    let mut i = 0;
    while i < NUM_OFFSET_SYMS {
        r[i] = OFFSETS[i].0 << 16 | OFFSETS[i].1;
        i += 1;
    }
    r
};

/// How a table's main part is sized.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum TableBits {
    /// Exactly this many bits.
    Fixed(u32),
    /// At most this many, fewer if the longest codeword is shorter: building a smaller table
    /// is faster, which matters for short blocks.
    AtMost(u32),
}

/// One entry: the symbol's static part plus `len` in both length fields.
fn entry(results: &[u32], sym: u16, len: u32) -> u32 {
    results[usize::from(sym)] + (len << 8) + len
}

/// Builds the decode table for the canonical Huffman code with codeword lengths `lens`
/// (each at most `max_len`) into `table`, and returns the main table's bits.
///
/// Returns `None` if the lengths are not a valid code: over-subscribed, or incomplete other
/// than in the two ways RFC 1951 allows (no codewords, or one of length 1). `sorted` is
/// scratch space for at least `lens.len()` symbols.
pub(super) fn build(
    table: &mut [u32],
    lens: &[u8],
    results: &[u32],
    bits: TableBits,
    mut max_len: u32,
    sorted: &mut [u16],
    pairs: Pairs,
) -> Option<u32> {
    // Neighbouring symbols often share a length, so counting them, and placing them in
    // order below, by one counter per length makes each step wait on the store of the last.
    // Instead, split the symbols into four contiguous quarters, each with its own counters,
    // and step through the four together.
    debug_assert!(
        lens.iter().all(|&l| u32::from(l) <= MAX_CODEWORD_LEN),
        "{lens:?}"
    );
    let q = lens.len() / 4;
    let quarters: [&[u8]; 4] = [
        &lens[..q],
        &lens[q..2 * q],
        &lens[2 * q..3 * q],
        &lens[3 * q..],
    ];
    let mut counts = [[0u32; 16]; 4];
    let [q0, q1, q2, q3] = quarters;
    for (((&a, &b), &c), &d) in q0.iter().zip(q1).zip(q2).zip(q3) {
        for (count, len) in counts.iter_mut().zip([a, b, c, d]) {
            count[usize::from(len & 15)] += 1;
        }
    }
    for &len in &quarters[3][q..] {
        counts[3][usize::from(len & 15)] += 1;
    }
    let mut len_counts = [0u32; MAX_CODEWORD_LEN as usize + 1];
    for (len, count) in len_counts.iter_mut().enumerate() {
        *count = counts.iter().map(|c| c[len]).sum();
    }
    while max_len > 1 && len_counts[max_len as usize] == 0 {
        max_len -= 1;
    }
    let table_bits = match bits {
        TableBits::Fixed(b) => b,
        TableBits::AtMost(b) => b.min(max_len),
    };

    // Sort the symbols by codeword length, then by value: canonical codeword order. Each
    // quarter writes its symbols of each length after the earlier quarters' ones.
    let mut next = [[0u32; 16]; 4];
    let mut start = 0;
    for len in 0..16 {
        for (n, c) in next.iter_mut().zip(&counts) {
            n[len] = start;
            start += c[len];
        }
    }
    for (j, (((&a, &b), &c), &d)) in q0.iter().zip(q1).zip(q2).zip(q3).enumerate() {
        for (k, (n, len)) in next.iter_mut().zip([a, b, c, d]).enumerate() {
            let slot = &mut n[usize::from(len & 15)];
            sorted[*slot as usize] = (k * q + j) as u16;
            *slot += 1;
        }
    }
    for (j, &len) in quarters[3].iter().enumerate().skip(q) {
        let slot = &mut next[3][usize::from(len & 15)];
        sorted[*slot as usize] = (3 * q + j) as u16;
        *slot += 1;
    }
    // The codespace the code uses, out of 2^max_len.
    let mut codespace_used = 0u32;
    for &count in &len_counts[1..=max_len as usize] {
        codespace_used = (codespace_used << 1) + count;
    }
    // Skip the unused symbols.
    let mut sorted = sorted[len_counts[0] as usize..lens.len()].iter().copied();

    if codespace_used > 1 << max_len {
        return None;
    }
    if codespace_used < 1 << max_len {
        // The allowed incomplete codes map to a complete one: one symbol (the used one, or
        // any) gets both 1-bit codewords.
        let sym = if codespace_used == 0 {
            0
        } else if codespace_used == 1 << (max_len - 1) && len_counts[1] == 1 {
            sorted.next()?
        } else {
            return None;
        };
        table[..1 << table_bits].fill(entry(results, sym, 1));
        return Some(table_bits);
    }

    // The code is complete. Fill the main table for the codewords of at most `table_bits`
    // bits, in codeword order. Codewords are bit-reversed, so one of length `len` repeats
    // every 2^len entries; rather than strided stores, treat the table as 2^len entries
    // while filling length `len`, and double it by copying as `len` grows. Within a length,
    // canonical codewords count up by one, so each one's reversed index is a lookup that
    // does not wait on the last.
    //
    // With `Pairs::Add`, literal pairs go in as well: each is written once, when the table is
    // as wide as its two codes together, and the doublings after copy it to every slot it
    // takes. A pair's slot starts with its first literal's code, so no other pair or
    // subtable pointer ever takes it: the code is prefix-free.
    let mut lits = [(0u32, 0u32); 256];
    let mut lit_start = [0usize; MAX_CODEWORD_LEN as usize + 2];
    let mut num_lits = 0;
    let add_pairs = pairs == Pairs::Add;
    let mut code = 0u32;
    let mut len = 1;
    let mut count = len_counts[1];
    while count == 0 {
        len += 1;
        count = len_counts[len as usize];
    }
    let mut cur_table_end = 1u32 << len;
    while len <= table_bits {
        for _ in 0..count {
            let sym = sorted.next()?;
            let rev = reversed(code, len);
            table[rev as usize] = entry(results, sym, len);
            if add_pairs && sym < 256 {
                lits[num_lits] = (rev, u32::from(sym));
                num_lits += 1;
            }
            code += 1;
        }
        lit_start[len as usize + 1] = num_lits;
        if add_pairs {
            add_pairs_of_span(table, &lits, &lit_start, len);
        }
        if code == 1 << len {
            // That was the last codeword, all ones.
            while len < table_bits {
                let end = cur_table_end as usize;
                table.copy_within(..end, end);
                cur_table_end <<= 1;
                len += 1;
                lit_start[len as usize + 1] = num_lits;
                if add_pairs {
                    add_pairs_of_span(table, &lits, &lit_start, len);
                }
            }
            return Some(table_bits);
        }
        loop {
            len += 1;
            code <<= 1;
            if len <= table_bits {
                let end = cur_table_end as usize;
                table.copy_within(..end, end);
                cur_table_end <<= 1;
            }
            count = len_counts[len as usize];
            if count != 0 {
                break;
            }
            // A width with no codewords of its own still gets its pairs before the next
            // doubling.
            lit_start[len as usize + 1] = num_lits;
            if add_pairs && len <= table_bits {
                add_pairs_of_span(table, &lits, &lit_start, len);
            }
        }
    }
    // Past the main table now, so `len` may exceed what `reversed` covers.
    let mut codeword = code.reverse_bits() >> (32 - len);

    // The longer codewords, in subtables after the main table.
    let prefix_mask = (1 << table_bits) - 1;
    cur_table_end = 1 << table_bits;
    let mut subtable_prefix = u32::MAX;
    let mut subtable_start = 0;
    loop {
        if codeword & prefix_mask != subtable_prefix {
            // A new subtable, big enough for the codewords of this length with this prefix,
            // and longer ones as needed to fill it (the code is complete, so they exist).
            subtable_prefix = codeword & prefix_mask;
            subtable_start = cur_table_end;
            let mut subtable_bits = len - table_bits;
            let mut used = count;
            while used < 1 << subtable_bits {
                subtable_bits += 1;
                used = (used << 1) + len_counts[(table_bits + subtable_bits) as usize];
            }
            cur_table_end = subtable_start + (1 << subtable_bits);
            table[subtable_prefix as usize] = subtable_start << 16
                | EXCEPTIONAL
                | SUBTABLE_POINTER
                | subtable_bits << 8
                | table_bits;
        }
        let e = entry(results, sorted.next()?, len - table_bits);
        let stride = 1 << (len - table_bits);
        let mut i = subtable_start + (codeword >> table_bits);
        while i < cur_table_end {
            table[i as usize] = e;
            i += stride;
        }
        if codeword == (1 << len) - 1 {
            return Some(table_bits);
        }
        codeword = next_codeword(codeword, len);
        count -= 1;
        while count == 0 {
            len += 1;
            count = len_counts[len as usize];
        }
    }
}

/// Whether [`build`] adds literal [`PAIR`] entries (only meaningful for the literal/length
/// table).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Pairs {
    None,
    Add,
}

/// Writes the pairs whose two literal codes take `span` bits together into the table, now
/// `2^span` entries wide: each once, at the codes' bit-reversed concatenation. `lits` holds
/// the literals' bit-reversed codes and symbols by code length, `lit_start[len]..
/// lit_start[len + 1]`, for every length below `span`.
fn add_pairs_of_span(table: &mut [u32], lits: &[(u32, u32)], lit_start: &[usize], span: u32) {
    for a in 1..span {
        let b = span - a;
        let firsts = &lits[lit_start[a as usize]..lit_start[a as usize + 1]];
        let seconds = &lits[lit_start[b as usize]..lit_start[b as usize + 1]];
        for &(rev1, sym1) in firsts {
            let entry = LITERAL | PAIR | sym1 << 16 | span;
            for &(rev2, sym2) in seconds {
                table[(rev1 | rev2 << a) as usize] = entry | sym2 << 8;
            }
        }
    }
}

/// `REVERSED[i]` is the 12-bit `i` with its bits in reverse order.
static REVERSED: [u16; 1 << LITLEN_TABLEBITS] = {
    let mut r = [0; 1 << LITLEN_TABLEBITS];
    let mut i = 0;
    while i < r.len() {
        r[i] = ((i as u32).reverse_bits() >> (32 - LITLEN_TABLEBITS)) as u16;
        i += 1;
    }
    r
};

/// The `len`-bit canonical codeword `code`, `len <= 12`, bit-reversed: DEFLATE's order, low
/// bit first.
fn reversed(code: u32, len: u32) -> u32 {
    debug_assert!(len <= LITLEN_TABLEBITS, "{len}-bit codeword");
    u32::from(REVERSED[(code as usize) & ((1 << LITLEN_TABLEBITS) - 1)]) >> (LITLEN_TABLEBITS - len)
}

/// The next codeword of `len` bits in canonical order, bit-reversed: increment the reversed
/// value, by setting its highest clear bit and clearing the bits above it.
fn next_codeword(codeword: u32, len: u32) -> u32 {
    let bit = 1 << (31 - (codeword ^ ((1 << len) - 1)).leading_zeros());
    (codeword & (bit - 1)) | bit
}
