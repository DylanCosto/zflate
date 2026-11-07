//! Deciding where to end a block, as libdeflate does: compare the mix of literals and
//! matches since the last check with the block's mix so far, and end the block when they
//! differ enough that new codes would pay for themselves.

use super::{
    block::Freqs,
    tables::{FIRST_LEN_SYM, NUM_LENGTH_SLOTS, NUM_LITERALS},
};

/// A block is at least this long, unless it is the last.
pub(super) const MIN_BLOCK_LENGTH: usize = 5000;

/// Literals by two of their high bits and their low bit; matches by whether they are
/// shorter than 9.
const NUM_LITERAL_OBSERVATION_TYPES: usize = 8;
pub(super) const NUM_OBSERVATION_TYPES: usize = NUM_LITERAL_OBSERVATION_TYPES + 2;
const OBSERVATIONS_PER_CHECK: u32 = 512;

pub(super) struct SplitStats {
    new_observations: [u32; NUM_OBSERVATION_TYPES],
    observations: [u32; NUM_OBSERVATION_TYPES],
    num_new_observations: u32,
    num_observations: u32,
}

impl SplitStats {
    pub(super) const fn new() -> Self {
        Self {
            new_observations: [0; NUM_OBSERVATION_TYPES],
            observations: [0; NUM_OBSERVATION_TYPES],
            num_new_observations: 0,
            num_observations: 0,
        }
    }

    pub(super) fn reset(&mut self) {
        *self = Self::new();
    }

    #[inline(always)]
    pub(super) fn observe_literal(&mut self, lit: u8) {
        self.new_observations[usize::from(((lit >> 5) & 0x6) | (lit & 1))] += 1;
        self.num_new_observations += 1;
    }

    #[inline(always)]
    pub(super) fn observe_match(&mut self, length: u32) {
        self.new_observations[NUM_LITERAL_OBSERVATION_TYPES + usize::from(length >= 9)] += 1;
        self.num_new_observations += 1;
    }

    /// Whether to end the block, for the greedy and lazy parsers, whose observations are
    /// exactly the symbols they choose: the new observations'
    /// types come from the symbol counts since the last check, so the parser only counts
    /// how many there are, `num_new`, in a register (reset here on a merge). Counting each
    /// type in memory, a read-modify-write per symbol, chains the symbols together through
    /// the store buffer; the decisions are the same.
    #[inline(always)]
    pub(super) fn should_end_block_from(
        &mut self,
        start: usize,
        pos: usize,
        end: usize,
        num_new: &mut u32,
        freqs: &Freqs,
    ) -> bool {
        *num_new >= OBSERVATIONS_PER_CHECK
            && pos - start >= MIN_BLOCK_LENGTH
            && end - pos >= MIN_BLOCK_LENGTH
            && self.check_from(num_new, freqs, (pos - start) as u32)
    }

    fn check_from(&mut self, num_new: &mut u32, freqs: &Freqs, block_len: u32) -> bool {
        let mut totals = [0; NUM_OBSERVATION_TYPES];
        for (lit, &freq) in freqs.litlen[..NUM_LITERALS].iter().enumerate() {
            totals[((lit >> 5) & 0x6) | (lit & 1)] += freq;
        }
        // Length slots 0 to 5 are lengths 3 to 8, one each.
        let lengths = &freqs.litlen[FIRST_LEN_SYM..FIRST_LEN_SYM + NUM_LENGTH_SLOTS];
        totals[NUM_LITERAL_OBSERVATION_TYPES] = lengths[..6].iter().sum();
        totals[NUM_LITERAL_OBSERVATION_TYPES + 1] = lengths[6..].iter().sum();
        for ((new, &total), &old) in self
            .new_observations
            .iter_mut()
            .zip(&totals)
            .zip(&self.observations)
        {
            *new = total - old;
        }
        self.num_new_observations = *num_new;
        let end_block = self.end_block_check(block_len);
        if !end_block {
            *num_new = 0;
        }
        end_block
    }

    /// Whether to check, every 512 observations, if the block that started at `start`
    /// should end before `pos`, with `end` the end of the input: once the block and what is
    /// left are both at least [`MIN_BLOCK_LENGTH`].
    #[inline(always)]
    pub(super) fn ready_to_check(&self, start: usize, pos: usize, end: usize) -> bool {
        self.num_new_observations >= OBSERVATIONS_PER_CHECK
            && pos - start >= MIN_BLOCK_LENGTH
            && end - pos >= MIN_BLOCK_LENGTH
    }

    /// Whether the new observations differ enough from the block's to end it before them;
    /// if not, merges them into the block's.
    pub(super) fn end_block_check(&mut self, block_len: u32) -> bool {
        if self.num_observations > 0 {
            // How far the new observations' distribution is from the block's, both scaled
            // to the same total.
            let total_delta: u32 = (0..NUM_OBSERVATION_TYPES)
                .map(|i| {
                    let expected = self.observations[i] * self.num_new_observations;
                    let actual = self.new_observations[i] * self.num_observations;
                    actual.abs_diff(expected)
                })
                .sum();
            let num_items = self.num_observations + self.num_new_observations;
            let mut cutoff = self.num_new_observations * 200 / 512 * self.num_observations;
            // Be more reluctant to end short blocks with few items, whose statistics are
            // noisier and whose codes cost relatively more.
            if block_len < 10000 && num_items < 8192 {
                cutoff += (u64::from(cutoff) * u64::from(8192 - num_items) / 8192) as u32;
            }
            if total_delta + (block_len / 4096) * self.num_observations >= cutoff {
                return true;
            }
        }
        self.merge();
        false
    }

    /// Adds the new observations to the block's.
    pub(super) fn merge(&mut self) {
        for (all, new) in self.observations.iter_mut().zip(&mut self.new_observations) {
            *all += *new;
            *new = 0;
        }
        self.num_observations += self.num_new_observations;
        self.num_new_observations = 0;
    }

    /// The block's observations, and how many there are.
    pub(super) fn observations(&self) -> (&[u32; NUM_OBSERVATION_TYPES], u32) {
        (&self.observations, self.num_observations)
    }

    /// Forgets the block's observations but keeps the new ones, which belong to the next
    /// block when a block is ended before them.
    pub(super) fn clear_old(&mut self) {
        self.observations = [0; NUM_OBSERVATION_TYPES];
        self.num_observations = 0;
    }
}
