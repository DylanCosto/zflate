//! Levels 10 to 12, as libdeflate's: find every useful match at every position with binary
//! trees, cache them, then choose the block's literals and matches by the cheapest path
//! through the block under a cost model, which is refined from the codes the path gives,
//! over several passes.

use super::{
    bitstream::BitWriter,
    block::{Body, OPTIMUM_LEN_MASK, OPTIMUM_OFFSET_SHIFT, OptimumNode, Sequence},
    default_costs::DEFAULT_COSTS,
    matchfinder::{BinaryTrees, LzMatch, TREES_REQUIRED_BYTES},
    parse::{Ctx, calculate_min_match_len, choose_min_match_len, clamp_lens, max_block_end},
    split::NUM_OBSERVATION_TYPES,
    tables::{
        END_OF_BLOCK, FIRST_LEN_SYM, LENGTH_EXTRA_BITS, LENGTH_SLOT, MAX_MATCH_LEN, MIN_MATCH_LEN,
        NUM_LITERALS, NUM_OFFSET_SLOTS, NUM_OFFSET_SYMS, OFFSET_EXTRA_BITS, OFFSET_SLOT_BASE,
    },
};

/// Costs are in sixteenths of a bit.
const BIT_COST: u32 = 16;
/// The costs, in bits, of symbols the last codes did not use.
const LITERAL_NOSTAT_BITS: u32 = 13;
const LENGTH_NOSTAT_BITS: u32 = 13;
const OFFSET_NOSTAT_BITS: u32 = 10;
const SOFT_MAX_BLOCK_LENGTH: usize = 300_000;
/// The match cache holds at most about this many entries before the block is ended.
const MATCH_CACHE_LENGTH: usize = SOFT_MAX_BLOCK_LENGTH * 5;
/// The most entries the match cache holds: [`MATCH_CACHE_LENGTH`] less one before the check that
/// ends the block, then a position's matches (up to 256, of increasing lengths), its count, and
/// the positions a nice match covers (up to 257).
const CACHE_CAPACITY: usize = MATCH_CACHE_LENGTH + 2 * MAX_MATCH_LEN as usize;
/// The most nodes a block takes: its length, at most [`SOFT_MAX_BLOCK_LENGTH`] and the positions
/// a nice match covers past it, and a match's length past its end.
const NODES_CAPACITY: usize = SOFT_MAX_BLOCK_LENGTH + 2 * MAX_MATCH_LEN as usize;
/// The cost to the end from a node past the block's end: too much for any path to use.
const COST_PAST_END: u32 = 0x8000_0000;

/// The slot of every offset, 1 to 32,768, for the path search's many lookups.
static OFFSET_SLOT_FULL: [u8; 32769] = offset_slots_full();

const fn offset_slots_full() -> [u8; 32769] {
    let mut table = [0; 32769];
    let mut slot = 0;
    while slot < NUM_OFFSET_SLOTS {
        let base = OFFSET_SLOT_BASE[slot] as usize;
        let mut offset = base;
        while offset < base + (1 << OFFSET_EXTRA_BITS[slot]) {
            table[offset] = slot as u8;
            offset += 1;
        }
        slot += 1;
    }
    table
}

/// How hard a level tries.
#[derive(Clone, Copy, Debug)]
pub(super) struct NearOptimalParams {
    pub(super) max_search_depth: u32,
    pub(super) nice_match_length: u32,
    /// Passes of the path search, each with costs from the codes of the pass before.
    pub(super) max_optim_passes: u32,
    /// Stop when a pass saves fewer bits than this.
    pub(super) min_improvement_to_continue: u32,
    /// Go back to an earlier pass's path when the last is worse by this many bits.
    pub(super) min_bits_to_use_nonfinal_path: u32,
    /// Blocks up to this long also get a path costed for the static codes.
    pub(super) max_len_to_optimize_static_block: usize,
}

/// The cost model: what each literal, length and offset slot costs.
#[derive(Clone)]
struct Costs {
    literal: [u32; NUM_LITERALS],
    length: [u32; MAX_MATCH_LEN as usize + 1],
    offset_slot: [u32; NUM_OFFSET_SYMS],
}

impl Costs {
    const fn zeroed() -> Self {
        Self {
            literal: [0; NUM_LITERALS],
            length: [0; MAX_MATCH_LEN as usize + 1],
            offset_slot: [0; NUM_OFFSET_SYMS],
        }
    }

    /// The costs of the codes with these lengths, unused symbols costing the `NOSTAT` bits.
    fn set_from_lens(&mut self, litlen: &[u8], offset: &[u8]) {
        let bits = |len: u8, nostat: u32| if len == 0 { nostat } else { u32::from(len) };
        for (cost, &len) in self.literal.iter_mut().zip(&litlen[..NUM_LITERALS]) {
            *cost = bits(len, LITERAL_NOSTAT_BITS) * BIT_COST;
        }
        for len in MIN_MATCH_LEN..=MAX_MATCH_LEN {
            let slot = usize::from(LENGTH_SLOT[len as usize]);
            let b = bits(litlen[FIRST_LEN_SYM + slot], LENGTH_NOSTAT_BITS)
                + u32::from(LENGTH_EXTRA_BITS[slot]);
            self.length[len as usize] = b * BIT_COST;
        }
        for slot in 0..NUM_OFFSET_SLOTS {
            let b = bits(offset[slot], OFFSET_NOSTAT_BITS) + u32::from(OFFSET_EXTRA_BITS[slot]);
            self.offset_slot[slot] = b * BIT_COST;
        }
    }
}

/// A length's starting cost, with its extra bits.
fn default_length_cost(len: u32, len_sym_cost: u32) -> u32 {
    len_sym_cost + u32::from(LENGTH_EXTRA_BITS[usize::from(LENGTH_SLOT[len as usize])]) * BIT_COST
}

/// An offset slot's starting cost: about 4.9 bits for the symbol, and its extra bits.
fn default_offset_slot_cost(slot: usize) -> u32 {
    4 * BIT_COST + (907 * BIT_COST) / 1000 + u32::from(OFFSET_EXTRA_BITS[slot]) * BIT_COST
}

pub(super) struct NearOptimal {
    trees: BinaryTrees,
    params: NearOptimalParams,
    /// The matches at each position of the block so far, each position's followed by their
    /// count and the position's byte.
    cache: Vec<LzMatch>,
    nodes: Vec<OptimumNode>,
    costs: Costs,
    costs_saved: Costs,
    /// The previous block's split observations, to judge how far its costs carry over.
    prev_observations: [u32; NUM_OBSERVATION_TYPES],
    prev_num_observations: u32,
    /// How often each length was the longest match, at the positions observed.
    new_match_len_freqs: [u32; MAX_MATCH_LEN as usize + 1],
    match_len_freqs: [u32; MAX_MATCH_LEN as usize + 1],
}

impl NearOptimal {
    /// The parser's tables, with room for the most matches and nodes a block takes, so that
    /// compressing allocates nothing; `None` if they cannot be allocated.
    pub(super) fn try_new(params: NearOptimalParams) -> Option<Self> {
        let mut cache = Vec::new();
        cache.try_reserve_exact(CACHE_CAPACITY).ok()?;
        let mut nodes = Vec::new();
        nodes.try_reserve_exact(NODES_CAPACITY).ok()?;
        Some(Self {
            trees: BinaryTrees::try_new()?,
            params,
            cache,
            nodes,
            costs: Costs::zeroed(),
            costs_saved: Costs::zeroed(),
            prev_observations: [0; NUM_OBSERVATION_TYPES],
            prev_num_observations: 0,
            new_match_len_freqs: [0; MAX_MATCH_LEN as usize + 1],
            match_len_freqs: [0; MAX_MATCH_LEN as usize + 1],
        })
    }

    fn init_stats(&mut self, ctx: &mut Ctx) {
        ctx.split.reset();
        self.new_match_len_freqs.fill(0);
        self.match_len_freqs.fill(0);
    }

    fn merge_stats(&mut self, ctx: &mut Ctx) {
        ctx.split.merge();
        for (all, new) in self
            .match_len_freqs
            .iter_mut()
            .zip(&mut self.new_match_len_freqs)
        {
            *all += *new;
            *new = 0;
        }
    }

    fn save_stats(&mut self, ctx: &Ctx) {
        let (observations, num) = ctx.split.observations();
        self.prev_observations = *observations;
        self.prev_num_observations = num;
    }

    fn clear_old_stats(&mut self, ctx: &mut Ctx) {
        ctx.split.clear_old();
        self.match_len_freqs.fill(0);
    }

    /// Compresses `data`: finds and caches the matches at each position, deciding where
    /// blocks end with the split statistics, and optimizes each block.
    #[inline(never)]
    pub(super) fn compress(&mut self, ctx: &mut Ctx, data: &[u8], w: &mut BitWriter<'_>) {
        let end = data.len();
        let (depth, nice) = (self.params.max_search_depth, self.params.nice_match_length);
        let mut pos = 0;
        let mut block_start = 0;
        let mut max_len = MAX_MATCH_LEN;
        let mut nice_len = nice.min(max_len);
        let mut prev_block_used_only_literals = false;
        self.trees.reset(data);
        self.cache.clear();
        self.init_stats(ctx);
        loop {
            let block_end_max = max_block_end(block_start, end, SOFT_MAX_BLOCK_LENGTH);
            let mut prev_end_block_check = None;
            let mut change_detected = false;
            let mut next_observation = pos;
            // After a block of only literals, observe only literals.
            let min_len = if prev_block_used_only_literals {
                MAX_MATCH_LEN + 1
            } else {
                calculate_min_match_len(&data[block_start..block_end_max], depth)
            };
            loop {
                let matches_start = self.cache.len();
                let mut best_len = 0;
                clamp_lens(&mut max_len, &mut nice_len, end - pos);
                if max_len >= TREES_REQUIRED_BYTES {
                    self.trees
                        .get_matches(data, pos, max_len, nice_len, depth, &mut self.cache);
                    if let Some(m) = self.cache[matches_start..].last() {
                        best_len = u32::from(m.length);
                    }
                }
                if pos >= next_observation {
                    if best_len >= min_len {
                        ctx.split.observe_match(best_len);
                        next_observation = pos + best_len as usize;
                        self.new_match_len_freqs[best_len as usize] += 1;
                    } else {
                        ctx.split.observe_literal(data[pos]);
                        next_observation = pos + 1;
                    }
                }
                self.cache.push(LzMatch {
                    length: (self.cache.len() - matches_start) as u16,
                    offset: u16::from(data[pos]),
                });
                pos += 1;
                // After a long match, only insert the positions it covers; a nice match is
                // almost surely the one to take.
                if best_len >= MIN_MATCH_LEN && best_len >= nice_len {
                    for _ in 1..best_len {
                        clamp_lens(&mut max_len, &mut nice_len, end - pos);
                        if max_len >= TREES_REQUIRED_BYTES {
                            self.trees.skip(data, pos, nice_len, depth);
                        }
                        self.cache.push(LzMatch {
                            length: 0,
                            offset: u16::from(data[pos]),
                        });
                        pos += 1;
                    }
                }
                if pos >= block_end_max || self.cache.len() >= MATCH_CACHE_LENGTH {
                    break;
                }
                if !ctx.split.ready_to_check(block_start, pos, end) {
                    continue;
                }
                if ctx.split.end_block_check((pos - block_start) as u32) {
                    change_detected = true;
                    break;
                }
                self.merge_stats(ctx);
                prev_end_block_check = Some(pos);
            }
            if let (true, Some(block_end)) = (change_detected, prev_end_block_check) {
                // The block changed character after the last check: end it there, and keep
                // the matches found since for the next block.
                let mut cache_end = self.cache.len();
                for _ in block_end..pos {
                    cache_end -= 1;
                    cache_end -= usize::from(self.cache[cache_end].length);
                }
                self.optimize_and_flush_block(
                    ctx,
                    w,
                    &data[block_start..block_end],
                    cache_end,
                    block_start == 0,
                    false,
                    &mut prev_block_used_only_literals,
                );
                self.cache.drain(..cache_end);
                self.save_stats(ctx);
                self.clear_old_stats(ctx);
                block_start = block_end;
            } else {
                self.merge_stats(ctx);
                let cache_end = self.cache.len();
                self.optimize_and_flush_block(
                    ctx,
                    w,
                    &data[block_start..pos],
                    cache_end,
                    block_start == 0,
                    pos == end,
                    &mut prev_block_used_only_literals,
                );
                self.cache.clear();
                self.save_stats(ctx);
                self.init_stats(ctx);
                block_start = pos;
            }
            // A change is only detected with a minimum block's length still to come.
            debug_assert!(block_start < end || pos == end);
            if pos == end || w.overflowed() {
                return;
            }
        }
    }

    /// Chooses the block's literals and matches and writes it: the cheapest path under
    /// costs refined over several passes, or only literals, or the path for the static
    /// codes, whichever costs least.
    #[expect(
        clippy::too_many_arguments,
        reason = "the block, where it is, and its result"
    )]
    fn optimize_and_flush_block(
        &mut self,
        ctx: &mut Ctx,
        w: &mut BitWriter<'_>,
        block: &[u8],
        cache_end: usize,
        is_first_block: bool,
        is_final_block: bool,
        used_only_literals: &mut bool,
    ) {
        let block_len = block.len();
        let coder = &mut ctx.coder;
        coder.choose_all_literals(block);
        let only_lits_cost = coder.true_cost();
        // Matches may reach past the block's end; the nodes there make that too costly.
        self.nodes.clear();
        self.nodes.resize(block_len + 1, OptimumNode::default());
        self.nodes.resize(
            block_len + MAX_MATCH_LEN as usize,
            OptimumNode {
                cost_to_end: COST_PAST_END,
                item: 0,
            },
        );

        let mut static_cost = u32::MAX;
        if block_len <= self.params.max_len_to_optimize_static_block {
            self.costs_saved.clone_from(&self.costs);
            let (litlen, offset) = coder.fixed_lens();
            self.costs.set_from_lens(litlen, offset);
            self.find_min_cost_path(ctx, block_len, cache_end);
            static_cost = self.nodes[0].cost_to_end / BIT_COST + 7;
            self.costs.clone_from(&self.costs_saved);
        }

        self.set_initial_costs(ctx, block, is_first_block);
        let mut best_true_cost = u32::MAX;
        let mut true_cost;
        let mut passes = self.params.max_optim_passes;
        loop {
            self.find_min_cost_path(ctx, block_len, cache_end);
            true_cost = ctx.coder.true_cost();
            if true_cost.saturating_add(self.params.min_improvement_to_continue) > best_true_cost {
                break;
            }
            best_true_cost = true_cost;
            self.costs_saved.clone_from(&self.costs);
            let (litlen, offset) = ctx.coder.lens();
            self.costs.set_from_lens(litlen, offset);
            passes -= 1;
            if passes == 0 {
                break;
            }
        }

        *used_only_literals = only_lits_cost < static_cost && only_lits_cost < best_true_cost;
        if *used_only_literals {
            ctx.coder.choose_all_literals(block);
            let (litlen, offset) = ctx.coder.lens();
            self.costs.set_from_lens(litlen, offset);
            let only_literals = [Sequence::last(block_len as u32)];
            ctx.coder
                .flush_block(w, block, Body::Sequences(&only_literals), is_final_block);
            return;
        }
        if static_cost < best_true_cost {
            // The path for the static codes; the block is still written with whichever
            // codes cost least for it.
            let (litlen, offset) = ctx.coder.fixed_lens();
            self.costs.set_from_lens(litlen, offset);
            self.find_min_cost_path(ctx, block_len, cache_end);
        } else if true_cost
            >= best_true_cost.saturating_add(self.params.min_bits_to_use_nonfinal_path)
        {
            // The last pass did worse: redo the best pass's path.
            self.costs.clone_from(&self.costs_saved);
            self.find_min_cost_path(ctx, block_len, cache_end);
            let (litlen, offset) = ctx.coder.lens();
            self.costs.set_from_lens(litlen, offset);
        }
        ctx.coder
            .flush_block(w, block, Body::Optimum(&self.nodes), is_final_block);
    }

    /// Finds the cheapest way to code the block under `self.costs`, from its end back to its
    /// start, using the cached matches (which end at `cache_end`); then counts the path's
    /// symbols and builds codes for them.
    fn find_min_cost_path(&mut self, ctx: &mut Ctx, block_len: usize, cache_end: usize) {
        let (costs, nodes, cache) = (&self.costs, &mut self.nodes, &self.cache);
        nodes[block_len].cost_to_end = 0;
        let mut c = cache_end;
        for i in (0..block_len).rev() {
            c -= 1;
            let header = cache[c];
            let num_matches = usize::from(header.length);
            let literal = u32::from(header.offset);
            let mut best = costs.literal[literal as usize] + nodes[i + 1].cost_to_end;
            let mut item = literal << OPTIMUM_OFFSET_SHIFT | 1;
            if num_matches > 0 {
                // The matches are in increasing order of length: each covers the lengths
                // from one past the previous one's.
                let mut len = MIN_MATCH_LEN as usize;
                for m in &cache[c - num_matches..c] {
                    let offset = u32::from(m.offset);
                    let offset_cost =
                        costs.offset_slot[usize::from(OFFSET_SLOT_FULL[offset as usize])];
                    while len <= usize::from(m.length) {
                        let cost = offset_cost + costs.length[len] + nodes[i + len].cost_to_end;
                        if cost < best {
                            best = cost;
                            item = len as u32 | offset << OPTIMUM_OFFSET_SHIFT;
                        }
                        len += 1;
                    }
                }
                c -= num_matches;
            }
            nodes[i] = OptimumNode {
                cost_to_end: best,
                item,
            };
        }
        // Count the path's symbols.
        let freqs = &mut ctx.coder.freqs;
        freqs.reset();
        let mut i = 0;
        while i < block_len {
            let item = nodes[i].item;
            let (len, offset) = (item & OPTIMUM_LEN_MASK, item >> OPTIMUM_OFFSET_SHIFT);
            if len == 1 {
                freqs.litlen[offset as usize] += 1;
            } else {
                freqs.litlen[FIRST_LEN_SYM + usize::from(LENGTH_SLOT[len as usize])] += 1;
                freqs.offset[usize::from(OFFSET_SLOT_FULL[offset as usize])] += 1;
            }
            i += len as usize;
        }
        freqs.litlen[END_OF_BLOCK] += 1;
        ctx.coder.make_codes();
    }

    /// The starting costs: libdeflate's defaults for the block's mix of literals and
    /// matches, or, after the first block, the last block's costs moved towards them as far
    /// as the block differs from the last.
    fn set_initial_costs(&mut self, ctx: &mut Ctx, block: &[u8], is_first_block: bool) {
        let (lit_cost, len_sym_cost) = self.default_litlen_costs(ctx, block);
        if is_first_block {
            self.set_default_costs(lit_cost, len_sym_cost);
        } else {
            self.adjust_costs(ctx, lit_cost, len_sym_cost);
        }
    }

    /// A literal's and a length symbol's default costs, from how many distinct literals the
    /// block uses and how often matches were found.
    fn default_litlen_costs(&self, ctx: &mut Ctx, block: &[u8]) -> (u32, u32) {
        let freqs = &mut ctx.coder.freqs.litlen[..NUM_LITERALS];
        freqs.fill(0);
        for &b in block {
            freqs[usize::from(b)] += 1;
        }
        let cutoff = block.len() as u32 >> 11;
        let num_used_literals = (freqs.iter().filter(|&&f| f > cutoff).count() as u32).max(1);
        // Estimate how many literals there will be from the matches found.
        let mut literal_freq = block.len() as i64;
        let mut match_freq = 0;
        let start = choose_min_match_len(num_used_literals, self.params.max_search_depth) as usize;
        for (len, &freq) in self.match_len_freqs.iter().enumerate().skip(start) {
            match_freq += i64::from(freq);
            literal_freq -= len as i64 * i64::from(freq);
        }
        let literal_freq = literal_freq.max(0);
        let i = if match_freq > literal_freq {
            2
        } else if match_freq * 4 > literal_freq {
            1
        } else {
            0
        };
        let costs = &DEFAULT_COSTS[i];
        (
            u32::from(costs.used_lits_to_lit_cost[num_used_literals as usize]),
            u32::from(costs.len_sym_cost),
        )
    }

    fn set_default_costs(&mut self, lit_cost: u32, len_sym_cost: u32) {
        self.costs.literal.fill(lit_cost);
        for len in MIN_MATCH_LEN..=MAX_MATCH_LEN {
            self.costs.length[len as usize] = default_length_cost(len, len_sym_cost);
        }
        for slot in 0..NUM_OFFSET_SLOTS {
            self.costs.offset_slot[slot] = default_offset_slot_cost(slot);
        }
    }

    /// Moves the costs towards the defaults, the more so the more this block's split
    /// observations differ from the last block's.
    fn adjust_costs(&mut self, ctx: &Ctx, lit_cost: u32, len_sym_cost: u32) {
        let (observations, num_observations) = ctx.split.observations();
        let total_delta: u64 = (0..NUM_OBSERVATION_TYPES)
            .map(|i| {
                let prev = u64::from(self.prev_observations[i]) * u64::from(num_observations);
                let cur = u64::from(observations[i]) * u64::from(self.prev_num_observations);
                prev.abs_diff(cur)
            })
            .sum();
        let cutoff =
            u64::from(self.prev_num_observations) * u64::from(num_observations) * 200 / 512;
        let change = if total_delta > 3 * cutoff {
            self.set_default_costs(lit_cost, len_sym_cost);
            return;
        } else if 4 * total_delta > 9 * cutoff {
            3
        } else if 2 * total_delta > 3 * cutoff {
            2
        } else if 2 * total_delta > cutoff {
            1
        } else {
            0
        };
        let adjust = |cost: &mut u32, default: u32| {
            *cost = match change {
                0 => (default + 3 * *cost) / 4,
                1 => (default + *cost) / 2,
                2 => (5 * default + 3 * *cost) / 8,
                _ => (3 * default + *cost) / 4,
            };
        };
        for cost in &mut self.costs.literal {
            adjust(cost, lit_cost);
        }
        for len in MIN_MATCH_LEN..=MAX_MATCH_LEN {
            adjust(
                &mut self.costs.length[len as usize],
                default_length_cost(len, len_sym_cost),
            );
        }
        for slot in 0..NUM_OFFSET_SLOTS {
            adjust(
                &mut self.costs.offset_slot[slot],
                default_offset_slot_cost(slot),
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::deflate::Deflater;

    /// Inputs that fill the match cache and the nodes most: several soft maximum blocks of
    /// random bytes, of one byte repeated (nice matches covering most positions), of a short
    /// period, and of text-like data, compressed at level 12, leave both as they were
    /// allocated.
    #[test]
    fn compressing_allocates_nothing() {
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let len = 3 * SOFT_MAX_BLOCK_LENGTH + 12_345;
        let inputs: [Vec<u8>; 4] = [
            (0..len).map(|_| next() as u8).collect(),
            vec![b'A'; len],
            (0..len).map(|i| b"ACGTTGCA"[i % 8]).collect(),
            (0..len)
                .map(|_| b"ACGT\t0123456789\n"[(next() % 16) as usize])
                .collect(),
        ];
        for level in [10, 11, 12] {
            let mut d = Deflater::try_new(level).expect("the tables fit");
            for input in &inputs {
                let mut out = vec![0; crate::deflate_compress_bound(input.len())];
                d.compress(input, &mut out).expect("fits");
                let super::super::Strategy::NearOptimal(n) = &d.strategy else {
                    unreachable!("levels 10 to 12 are near-optimal")
                };
                assert_eq!(n.cache.capacity(), CACHE_CAPACITY, "level {level}");
                assert_eq!(n.nodes.capacity(), NODES_CAPACITY, "level {level}");
            }
        }
    }
}
