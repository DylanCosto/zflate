//! Choosing the literals and matches, as libdeflate's levels 1 to 9 do: greedily (take the
//! longest match at each position), or lazily (first check whether a longer match starts
//! one or two bytes on).

use crate::alloc::try_boxed_array;

use super::{
    bitstream::BitWriter,
    block::{BlockCoder, Body, Freqs, Sequence},
    matchfinder::{HashChains, HashTable, Match, SMALL_INPUT_MAX, TABLE_REQUIRED_BYTES},
    split::{MIN_BLOCK_LENGTH, SplitStats},
    tables::{MAX_MATCH_LEN, MIN_MATCH_LEN, NUM_LITERALS, offset_slot},
};

/// The longest a block may get before it is ended regardless of its statistics, and the
/// most sequences it may hold; level 1's are smaller.
const SOFT_MAX_BLOCK_LENGTH: usize = 300_000;
pub(super) const SEQ_STORE_LENGTH: usize = 50_000;
const FAST_SOFT_MAX_BLOCK_LENGTH: usize = 65_535;
pub(super) const FAST_SEQ_STORE_LENGTH: usize = 8192;

/// How hard the hash-chain levels search.
#[derive(Clone, Copy, Debug)]
pub(super) struct SearchParams {
    pub(super) max_search_depth: u32,
    pub(super) nice_match_length: u32,
    /// The lazy parsers' depth for the search one byte on.
    pub(super) lookahead_depth: u32,
}

/// What the parsers keep between blocks and calls: the block coder and the split
/// statistics.
pub(super) struct Ctx {
    pub(super) coder: BlockCoder,
    pub(super) split: SplitStats,
}

/// Room for a block's sequences: at least the most any level keeps and the last one's
/// literals, rounded up to a power of two, so that an index masked to it needs no check.
pub(super) const SEQS_LEN: usize = (SEQ_STORE_LENGTH + 1).next_power_of_two();
pub(super) type Seqs = [Sequence; SEQS_LEN];

/// Allocates room for a block's sequences, or returns `None` if it cannot.
pub(super) fn try_new_seqs() -> Option<Box<Seqs>> {
    try_boxed_array(Sequence::default())
}

impl Ctx {
    pub(super) fn new() -> Self {
        Self {
            coder: BlockCoder::new(),
            split: SplitStats::new(),
        }
    }

    /// Starts a block, resetting its symbol counts.
    fn begin_block<'a>(&'a mut self, seqs: &'a mut Seqs) -> Block<'a> {
        self.coder.freqs.reset();
        Block {
            ctx: self,
            seqs,
            num_seqs: 0,
            litrunlen: 0,
            new_observations: 0,
        }
    }

    /// Ends the block `data`, of `num_seqs` sequences and the literals after them.
    fn end_block(
        &mut self,
        w: &mut BitWriter<'_>,
        data: &[u8],
        seqs: &mut Seqs,
        (num_seqs, litrunlen): (usize, u32),
        is_final: bool,
    ) {
        seqs[num_seqs] = Sequence::last(litrunlen);
        self.coder
            .finish_block(w, data, Body::Sequences(&seqs[..=num_seqs]), is_final);
    }
}

/// The block being parsed. A local, so that its counts stay in registers; it holds two
/// pointers, which the parsers, short of registers, can afford.
struct Block<'a> {
    ctx: &'a mut Ctx,
    seqs: &'a mut Seqs,
    num_seqs: usize,
    /// Literals since the last match.
    litrunlen: u32,
    /// Symbols since the last block-split check ([`SplitStats::should_end_block_from`]).
    new_observations: u32,
}

impl Block<'_> {
    #[inline(always)]
    fn literal<const SPLIT: bool>(&mut self, lit: u8) {
        self.ctx.coder.freqs.add_literal(lit);
        if SPLIT {
            self.new_observations += 1;
        }
        self.litrunlen += 1;
    }

    #[inline(always)]
    fn add_match<const SPLIT: bool>(&mut self, m: Match) {
        debug_assert!(
            m.len >= MIN_MATCH_LEN && m.len <= MAX_MATCH_LEN && (1..=32768).contains(&m.offset)
        );
        let slot = offset_slot(m.offset);
        self.ctx.coder.freqs.add_match(m.len, slot);
        if SPLIT {
            self.new_observations += 1;
        }
        // The parsers end a block before `SEQ_STORE_LENGTH` sequences, so the mask never
        // changes the index.
        self.seqs[self.num_seqs & (SEQS_LEN - 1)] =
            Sequence::new(self.litrunlen, m.len, m.offset, slot);
        self.num_seqs += 1;
        self.litrunlen = 0;
    }

    /// Whether to end the block that started at `start` before `pos`.
    #[inline(always)]
    fn should_end(&mut self, start: usize, pos: usize, end: usize) -> bool {
        let ctx = &mut *self.ctx;
        ctx.split.should_end_block_from(
            start,
            pos,
            end,
            &mut self.new_observations,
            &ctx.coder.freqs,
        )
    }

    /// The sequences and the literals after them, for [`Ctx::end_block`].
    fn finish(self) -> (usize, u32) {
        (self.num_seqs, self.litrunlen)
    }
}

/// `data[pos]`, or 0 past the end, which no caller reaches. A read that cannot panic keeps
/// the parsers' loops free of panic paths, which cost them more than the compare.
#[inline(always)]
fn byte_at(data: &[u8], pos: usize) -> u8 {
    data.get(pos).copied().unwrap_or(0)
}

/// Where a block starting at `start` must end at the latest: `soft_max` on, unless less
/// than a minimum block would be left after that.
pub(super) fn max_block_end(start: usize, end: usize, soft_max: usize) -> usize {
    if end - start < soft_max + MIN_BLOCK_LENGTH {
        end
    } else {
        start + soft_max
    }
}

/// Lowers the longest and nice match lengths near the end of the input.
#[inline(always)]
pub(super) fn clamp_lens(max_len: &mut u32, nice_len: &mut u32, remaining: usize) {
    if remaining < MAX_MATCH_LEN as usize {
        crate::cold_path();
        *max_len = remaining as u32;
        *nice_len = (*nice_len).min(*max_len);
    }
}

/// The shortest match worth taking, given how many distinct literals the data uses (the
/// table is libdeflate's, from its source): with
/// few, literals are cheap and short matches do not pay. Shallow searches find fewer long
/// matches, so they lower it.
pub(super) fn choose_min_match_len(num_used_literals: u32, max_search_depth: u32) -> u32 {
    #[rustfmt::skip]
    const MIN_LENS: [u8; 81] = [
        9, 9, 9, 9, 9, 9, 8, 8, 7, 7, 6, 6, 6, 6, 6, 6,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        3,
    ];
    let Some(&min_len) = MIN_LENS.get(num_used_literals as usize) else {
        return MIN_MATCH_LEN;
    };
    let min_len = u32::from(min_len);
    match max_search_depth {
        0..5 => min_len.min(4),
        5..10 => min_len.min(5),
        10..16 => min_len.min(7),
        _ => min_len,
    }
}

/// [`choose_min_match_len`] from the distinct bytes in the block's first 4 KiB.
pub(super) fn calculate_min_match_len(data: &[u8], max_search_depth: u32) -> u32 {
    if data.len() < 512 {
        return MIN_MATCH_LEN;
    }
    let mut used = [false; 256];
    for &b in &data[..data.len().min(4096)] {
        used[usize::from(b)] = true;
    }
    let num_used = used.iter().filter(|&&u| u).count() as u32;
    choose_min_match_len(num_used, max_search_depth)
}

/// [`choose_min_match_len`] from the literals chosen so far, counting those used more than
/// once in 1024 literals.
fn recalculate_min_match_len(freqs: &Freqs, max_search_depth: u32) -> u32 {
    let literals = &freqs.litlen[..NUM_LITERALS];
    let cutoff = literals.iter().sum::<u32>() >> 10;
    let num_used = literals.iter().filter(|&&f| f > cutoff).count() as u32;
    choose_min_match_len(num_used, max_search_depth)
}

/// The integer log2 of a nonzero offset.
#[inline(always)]
fn bsr32(x: u32) -> i32 {
    31 - x.leading_zeros() as i32
}

/// Level 1: the longest match in a 2-entry hash bucket, taken greedily; small blocks, no
/// split statistics.
pub(super) fn compress_fastest(
    table: &mut HashTable,
    ctx: &mut Ctx,
    seqs: &mut Seqs,
    nice_match_length: u32,
    data: &[u8],
    w: &mut BitWriter<'_>,
) {
    if data.len() <= SMALL_INPUT_MAX {
        compress_fastest_inner::<true>(table, ctx, seqs, nice_match_length, data, w);
    } else {
        compress_fastest_inner::<false>(table, ctx, seqs, nice_match_length, data, w);
    }
}

#[inline(never)]
fn compress_fastest_inner<const SMALL: bool>(
    table: &mut HashTable,
    ctx: &mut Ctx,
    seqs: &mut Seqs,
    nice_match_length: u32,
    data: &[u8],
    w: &mut BitWriter<'_>,
) {
    let end = data.len();
    let mut pos = 0;
    let mut max_len = MAX_MATCH_LEN;
    let mut nice_len = nice_match_length.min(max_len);
    let mut next_key = table.reset::<SMALL>(data);
    loop {
        let start = pos;
        let block_end = max_block_end(pos, end, FAST_SOFT_MAX_BLOCK_LENGTH);
        let mut b = ctx.begin_block(seqs);
        while pos < block_end && b.num_seqs < FAST_SEQ_STORE_LENGTH {
            let remaining = end - pos;
            if remaining < MAX_MATCH_LEN as usize {
                crate::cold_path();
                max_len = remaining as u32;
                if remaining < TABLE_REQUIRED_BYTES {
                    for &lit in &data[pos..] {
                        b.literal::<false>(lit);
                    }
                    pos = end;
                    break;
                }
                nice_len = nice_len.min(max_len);
            }
            // SAFETY: the tail check guarantees max_len bytes and at least five.
            // The caller uses SMALL only for inputs with 16-bit positions.
            let (len, offset) = unsafe {
                table.longest_match::<SMALL>(&mut next_key, data, pos, max_len, nice_len)
            };
            if len > 0 {
                b.add_match::<false>(Match { len, offset });
                table.skip::<SMALL>(&mut next_key, data, pos + 1, len as usize - 1, offset);
                pos += len as usize;
            } else {
                b.literal::<false>(byte_at(data, pos));
                pos += 1;
            }
        }
        let counts = b.finish();
        ctx.end_block(w, &data[start..pos], seqs, counts, pos == end);
        if pos == end || w.overflowed() {
            return;
        }
    }
}

/// Levels 2 to 4: the longest match on the hash chains, taken greedily.
#[inline(never)]
pub(super) fn compress_greedy<const SMALL: bool>(
    chains: &mut HashChains,
    ctx: &mut Ctx,
    seqs: &mut Seqs,
    params: SearchParams,
    data: &[u8],
    w: &mut BitWriter<'_>,
) {
    let end = data.len();
    let mut pos = 0;
    let mut max_len = MAX_MATCH_LEN;
    let mut nice_len = params.nice_match_length.min(max_len);
    let depth = params.max_search_depth;
    let mut hashes = chains.reset::<SMALL, true>(data);
    loop {
        let start = pos;
        let block_end = max_block_end(pos, end, SOFT_MAX_BLOCK_LENGTH);
        ctx.split.reset();
        let mut b = ctx.begin_block(seqs);
        let min_len = calculate_min_match_len(&data[pos..block_end], depth);
        while pos < block_end && b.num_seqs < SEQ_STORE_LENGTH && !b.should_end(start, pos, end) {
            clamp_lens(&mut max_len, &mut nice_len, end - pos);
            let m = chains.longest_match::<SMALL, true>(
                &mut hashes,
                data,
                pos,
                min_len - 1,
                max_len,
                nice_len,
                depth,
            );
            // A 3-byte match far back costs about as much as its literals.
            if m.len >= min_len && (m.len > MIN_MATCH_LEN || m.offset <= 4096) {
                b.add_match::<true>(m);
                chains.skip::<SMALL, true>(&mut hashes, data, pos + 1, m.len as usize - 1);
                pos += m.len as usize;
            } else {
                b.literal::<true>(byte_at(data, pos));
                pos += 1;
            }
        }
        let counts = b.finish();
        ctx.end_block(w, &data[start..pos], seqs, counts, pos == end);
        if pos == end || w.overflowed() {
            return;
        }
    }
}

/// Levels 5 to 9: before taking a match, look for a better one starting a byte later (and,
/// with `LAZY2`, two bytes later), and take that instead, after the skipped literals, if
/// it is longer by enough to pay for them and for any larger offset. Without `FIND3`, matches are
/// at least 4 bytes long.
#[inline(never)]
pub(super) fn compress_lazy<const LAZY2: bool, const SMALL: bool, const FIND3: bool>(
    chains: &mut HashChains,
    ctx: &mut Ctx,
    seqs: &mut Seqs,
    params: SearchParams,
    data: &[u8],
    w: &mut BitWriter<'_>,
) {
    let end = data.len();
    let mut pos = 0;
    let mut max_len = MAX_MATCH_LEN;
    let mut nice_len = params.nice_match_length.min(max_len);
    let depth = params.max_search_depth;
    // `skip_match` leaves out inserts that no search `depth` deep can reach.
    debug_assert!(params.lookahead_depth <= depth);
    let shortest = if FIND3 {
        MIN_MATCH_LEN
    } else {
        MIN_MATCH_LEN + 1
    };
    let mut hashes = chains.reset::<SMALL, FIND3>(data);
    loop {
        let start = pos;
        let block_end = max_block_end(pos, end, SOFT_MAX_BLOCK_LENGTH);
        let mut next_recalc_min_len = pos + (end - pos).min(10000);
        ctx.split.reset();
        let mut b = ctx.begin_block(seqs);
        let mut min_len = calculate_min_match_len(&data[pos..block_end], depth).max(shortest);
        while pos < block_end && b.num_seqs < SEQ_STORE_LENGTH && !b.should_end(start, pos, end) {
            // Revisit the minimum length from the literals chosen, each time the block
            // doubles.
            if pos >= next_recalc_min_len {
                min_len = recalculate_min_match_len(&b.ctx.coder.freqs, depth).max(shortest);
                next_recalc_min_len += (end - next_recalc_min_len).min(pos - start);
            }
            clamp_lens(&mut max_len, &mut nice_len, end - pos);
            let mut cur = chains.longest_match::<SMALL, FIND3>(
                &mut hashes,
                data,
                pos,
                min_len - 1,
                max_len,
                nice_len,
                depth,
            );
            if cur.len < min_len || (FIND3 && cur.len == MIN_MATCH_LEN && cur.offset > 8192) {
                b.literal::<true>(byte_at(data, pos));
                pos += 1;
                continue;
            }
            // `pos` is now one past the start of the match `cur`.
            pos += 1;
            loop {
                if cur.len >= nice_len {
                    b.add_match::<true>(cur);
                    chains.skip_match::<SMALL, FIND3>(
                        &mut hashes,
                        data,
                        pos,
                        cur.len as usize - 1,
                        cur.offset,
                        depth,
                    );
                    pos += cur.len as usize - 1;
                    break;
                }
                clamp_lens(&mut max_len, &mut nice_len, end - pos);
                let next = chains.longest_match::<SMALL, FIND3>(
                    &mut hashes,
                    data,
                    pos,
                    cur.len - 1,
                    max_len,
                    nice_len,
                    params.lookahead_depth,
                );
                pos += 1;
                if next.len >= cur.len
                    && 4 * (next.len - cur.len) as i32 + (bsr32(cur.offset) - bsr32(next.offset))
                        > 2
                {
                    b.literal::<true>(byte_at(data, pos - 2));
                    cur = next;
                    continue;
                }
                if LAZY2 {
                    clamp_lens(&mut max_len, &mut nice_len, end - pos);
                    let next = chains.longest_match::<SMALL, FIND3>(
                        &mut hashes,
                        data,
                        pos,
                        cur.len - 1,
                        max_len,
                        nice_len,
                        depth >> 2,
                    );
                    pos += 1;
                    if next.len >= cur.len
                        && 4 * (next.len - cur.len) as i32
                            + (bsr32(cur.offset) - bsr32(next.offset))
                            > 6
                    {
                        b.literal::<true>(byte_at(data, pos - 3));
                        b.literal::<true>(byte_at(data, pos - 2));
                        cur = next;
                        continue;
                    }
                    b.add_match::<true>(cur);
                    if cur.len > 3 {
                        chains.skip::<SMALL, FIND3>(&mut hashes, data, pos, cur.len as usize - 3);
                        pos += cur.len as usize - 3;
                    }
                } else {
                    b.add_match::<true>(cur);
                    chains.skip_match::<SMALL, FIND3>(
                        &mut hashes,
                        data,
                        pos,
                        cur.len as usize - 2,
                        cur.offset,
                        depth,
                    );
                    pos += cur.len as usize - 2;
                }
                break;
            }
        }
        let counts = b.finish();
        ctx.end_block(w, &data[start..pos], seqs, counts, pos == end);
        if pos == end || w.overflowed() {
            return;
        }
    }
}
