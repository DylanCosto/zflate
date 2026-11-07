//! DEFLATE's alphabets (RFC 1951 3.2.5) as the compressor needs them: slots for lengths and
//! offsets, their extra bits, and the static codes' lengths.

pub(super) const MIN_MATCH_LEN: u32 = 3;
pub(super) const MAX_MATCH_LEN: u32 = 258;
pub(super) const NUM_LITLEN_SYMS: usize = 288;
pub(super) const NUM_OFFSET_SYMS: usize = 32;
pub(super) const NUM_PRECODE_SYMS: usize = 19;
pub(super) const NUM_LITERALS: usize = 256;
pub(super) const END_OF_BLOCK: usize = 256;
pub(super) const FIRST_LEN_SYM: usize = 257;
pub(super) const NUM_LENGTH_SLOTS: usize = 29;
pub(super) const NUM_OFFSET_SLOTS: usize = 30;
/// The offset slot tables' length: DEFLATE's 30 slots, padded to 32 as the offset code is,
/// so that a slot masked to 5 bits needs no bounds check.
pub(super) const OFFSET_SLOT_TABLE_LEN: usize = 32;

/// The longest codewords each code may have. The literal/length code stops at 14 bits rather
/// than DEFLATE's 15, as libdeflate's does, so that a length's codeword and its extra bits
/// fit in 19 bits; that costs almost nothing.
pub(super) const MAX_LITLEN_CODEWORD_LEN: u32 = 14;
pub(super) const MAX_OFFSET_CODEWORD_LEN: u32 = 15;
pub(super) const MAX_PRE_CODEWORD_LEN: u32 = 7;

pub(super) const LENGTH_SLOT_BASE: [u32; NUM_LENGTH_SLOTS] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131,
    163, 195, 227, 258,
];
pub(super) const LENGTH_EXTRA_BITS: [u8; NUM_LENGTH_SLOTS] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];
pub(super) const OFFSET_SLOT_BASE: [u32; OFFSET_SLOT_TABLE_LEN] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
    2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577, 32769, 49153,
];
pub(super) const OFFSET_EXTRA_BITS: [u8; OFFSET_SLOT_TABLE_LEN] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13, 14, 14,
];
/// The order in which a dynamic block's header gives the precode's lengths.
pub(super) const PRECODE_LENS_PERMUTATION: [u8; NUM_PRECODE_SYMS] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];
pub(super) const PRECODE_EXTRA_BITS: [u8; NUM_PRECODE_SYMS] =
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 7];

/// Tables indexed by a match length have this many entries, padded past 258, so that a
/// length masked to 9 bits needs no bounds check.
pub(super) const LENGTH_TABLE_LEN: usize = 512;

/// The slot of each match length, 3 to 258.
pub(super) static LENGTH_SLOT: [u8; LENGTH_TABLE_LEN] = length_slots();

/// `LENGTH_SLOT[len]`, without a bounds check.
#[inline(always)]
pub(super) fn length_slot(len: u32) -> usize {
    usize::from(LENGTH_SLOT[len as usize & (LENGTH_TABLE_LEN - 1)])
}

const fn length_slots() -> [u8; LENGTH_TABLE_LEN] {
    let mut table = [0; LENGTH_TABLE_LEN];
    let mut slot = 0;
    let mut len = MIN_MATCH_LEN;
    while len <= MAX_MATCH_LEN {
        while slot + 1 < NUM_LENGTH_SLOTS && LENGTH_SLOT_BASE[slot + 1] <= len {
            slot += 1;
        }
        table[len as usize] = slot as u8;
        len += 1;
    }
    table
}

/// The slot of offsets 1 to 256, at `offset - 1`; [`offset_slot`] derives the rest from it.
static OFFSET_SLOT_SMALL: [u8; 256] = offset_slots_small();

const fn offset_slots_small() -> [u8; 256] {
    let mut table = [0; 256];
    let mut slot = 0;
    let mut offset = 1;
    while offset <= 256 {
        while slot + 1 < NUM_OFFSET_SLOTS && OFFSET_SLOT_BASE[slot + 1] <= offset {
            slot += 1;
        }
        table[offset as usize - 1] = slot as u8;
        offset += 1;
    }
    table
}

/// The slot of a match offset, 1 to 32,768. Slots 16 to 29 are each 128 times as wide as
/// slots 2 to 15, so an offset past 256 is looked up at `(offset - 1) >> 7` and moved up 14
/// slots.
#[inline(always)]
pub(super) fn offset_slot(offset: u32) -> usize {
    debug_assert!((1..=32768).contains(&offset));
    let n = if offset <= 256 { 0 } else { 7 };
    (usize::from(OFFSET_SLOT_SMALL[((offset - 1) >> n) as usize & 0xff]) + (n << 1) as usize) & 31
}

/// The static literal/length code's lengths (RFC 1951 3.2.6).
pub(super) const fn static_litlen_len(sym: usize) -> u8 {
    match sym {
        0..144 => 8,
        144..256 => 9,
        256..280 => 7,
        _ => 8,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slots_cover_their_ranges() {
        for len in MIN_MATCH_LEN..=MAX_MATCH_LEN {
            let slot = usize::from(LENGTH_SLOT[len as usize]);
            let base = LENGTH_SLOT_BASE[slot];
            assert!(
                base <= len && len - base < 1 << LENGTH_EXTRA_BITS[slot],
                "length {len}"
            );
            // 258 has its own slot, though slot 27 could also code it.
            assert_eq!(len == 258, slot == 28, "length {len}");
        }
        for offset in 1..=32768 {
            let slot = offset_slot(offset);
            let base = OFFSET_SLOT_BASE[slot];
            assert!(
                base <= offset && offset - base < 1 << OFFSET_EXTRA_BITS[slot],
                "offset {offset}"
            );
        }
    }
}
