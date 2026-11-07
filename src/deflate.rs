//! The DEFLATE compressor: libdeflate's design, ported. Each level picks a
//! match finder and a parser; the blocks they produce are ended by [`block`], which writes
//! each as whichever of a dynamic, static, or stored block is smallest.

mod bitstream;
mod block;
mod default_costs;
mod huffman;
mod matchfinder;
mod near_optimal;
mod parse;
mod split;
mod tables;

use crate::alloc::try_box;

use self::{
    bitstream::BitWriter,
    matchfinder::{HashChains, HashTable, SMALL_INPUT_MAX},
    near_optimal::{NearOptimal, NearOptimalParams},
    parse::{Ctx, SearchParams, Seqs, try_new_seqs},
    split::MIN_BLOCK_LENGTH,
};

/// Level 5's search depth, look-ahead depth and nice length.
const LEVEL5: (u32, u32, u32) = (12, 3, 48);

/// How a level compresses.
enum Strategy {
    /// Level 0: stored blocks only.
    Stored,
    /// Level 1.
    Fastest(HashTable),
    /// Levels 2 to 4.
    Greedy(HashChains),
    /// Level 5: lazy, but taking no 3-byte matches, and so keeping no 3-byte hash table. On
    /// BAM and BCF, libdeflate's 3-byte matches (within 8 KiB, when the input's bytes are
    /// varied) cost more than they save; text gets a longer minimum anyway.
    Lazy4(HashChains),
    /// Levels 6 and 7.
    Lazy(HashChains),
    /// Levels 8 and 9.
    Lazy2(HashChains),

    /// Levels 10 to 12.
    NearOptimal(Box<NearOptimal>),
}

/// A compressor for one level, with its match finder and block state, reused across inputs.
pub(crate) struct Deflater {
    strategy: Strategy,
    params: SearchParams,
    /// Inputs this short are stored: compressing them would rarely pay. Higher levels try
    /// harder.
    max_passthrough_size: usize,
    ctx: Ctx,
    /// Room for a block's sequences, for the levels that keep them.
    seqs: Option<Box<Seqs>>,
}

impl Deflater {
    /// A compressor at `level`, 0 to 12, or `None` if its tables cannot be allocated.
    pub(crate) fn try_new(level: u8) -> Option<Self> {
        let near_optimal =
            |max_search_depth, nice_match_length, passes, min_improvement, max_static| {
                let params = NearOptimalParams {
                    max_search_depth,
                    nice_match_length,
                    max_optim_passes: passes,
                    min_improvement_to_continue: min_improvement,
                    min_bits_to_use_nonfinal_path: min_improvement,
                    max_len_to_optimize_static_block: max_static,
                };
                Some(Strategy::NearOptimal(try_box(NearOptimal::try_new(
                    params,
                )?)?))
            };
        let (strategy, max_search_depth, nice_match_length) = match level {
            0 => (Strategy::Stored, 0, 0),
            1 => (Strategy::Fastest(HashTable::try_new()?), 0, 32),
            2 => (Strategy::Greedy(HashChains::try_new()?), 6, 10),
            3 => (Strategy::Greedy(HashChains::try_new()?), 12, 14),
            4 => (Strategy::Greedy(HashChains::try_new()?), 16, 30),
            5 => (Strategy::Lazy4(HashChains::try_new()?), LEVEL5.0, LEVEL5.2),
            6 => (Strategy::Lazy(HashChains::try_new()?), 35, 65),
            7 => (Strategy::Lazy(HashChains::try_new()?), 100, 130),
            8 => (Strategy::Lazy2(HashChains::try_new()?), 300, 258),
            9 => (Strategy::Lazy2(HashChains::try_new()?), 600, 258),
            10 => (near_optimal(35, 75, 2, 32, 0)?, 35, 75),
            11 => (near_optimal(100, 150, 4, 16, 1000)?, 100, 150),
            _ => (near_optimal(300, 258, 10, 1, 10_000)?, 300, 258),
        };
        let seqs = match strategy {
            Strategy::Stored | Strategy::NearOptimal(_) => None,
            _ => Some(try_new_seqs()?),
        };
        Some(Self {
            strategy,
            params: SearchParams {
                max_search_depth,
                nice_match_length,
                // libdeflate's look-ahead searches half as deep, but for level 5.
                lookahead_depth: if level == 5 {
                    LEVEL5.1
                } else {
                    max_search_depth >> 1
                },
            },
            max_passthrough_size: if level == 0 {
                usize::MAX
            } else {
                55 - 4 * usize::from(level.min(12))
            },
            ctx: Ctx::new(),
            seqs,
        })
    }

    /// Compresses `data` into `out` as raw DEFLATE, returning the length written, or `None`
    /// if it does not fit.
    pub(crate) fn compress(&mut self, data: &[u8], out: &mut [u8]) -> Option<usize> {
        if data.len() <= self.max_passthrough_size {
            return compress_stored(data, out);
        }
        let mut w = BitWriter::new(out);
        let (ctx, params) = (&mut self.ctx, self.params);
        let seqs = self.seqs.as_deref_mut();
        let seqs = || seqs.expect("the parsing levels keep sequences");
        match &mut self.strategy {
            Strategy::Stored => unreachable!("level 0 stores every input"),
            Strategy::Fastest(table) => {
                parse::compress_fastest(table, ctx, seqs(), params.nice_match_length, data, &mut w);
            }
            // Inputs up to 64 KiB, BGZF's, keep positions without a window to move.
            Strategy::Greedy(chains) if data.len() <= SMALL_INPUT_MAX => {
                parse::compress_greedy::<true>(chains, ctx, seqs(), params, data, &mut w);
            }
            Strategy::Greedy(chains) => {
                parse::compress_greedy::<false>(chains, ctx, seqs(), params, data, &mut w);
            }
            Strategy::Lazy4(chains) if data.len() <= SMALL_INPUT_MAX => {
                parse::compress_lazy::<false, true, false>(
                    chains,
                    ctx,
                    seqs(),
                    params,
                    data,
                    &mut w,
                );
            }
            Strategy::Lazy4(chains) => {
                parse::compress_lazy::<false, false, false>(
                    chains,
                    ctx,
                    seqs(),
                    params,
                    data,
                    &mut w,
                );
            }
            Strategy::Lazy(chains) if data.len() <= SMALL_INPUT_MAX => {
                parse::compress_lazy::<false, true, true>(
                    chains,
                    ctx,
                    seqs(),
                    params,
                    data,
                    &mut w,
                );
            }
            Strategy::Lazy(chains) => {
                parse::compress_lazy::<false, false, true>(
                    chains,
                    ctx,
                    seqs(),
                    params,
                    data,
                    &mut w,
                );
            }
            Strategy::Lazy2(chains) if data.len() <= SMALL_INPUT_MAX => {
                parse::compress_lazy::<true, true, true>(chains, ctx, seqs(), params, data, &mut w);
            }
            Strategy::Lazy2(chains) => {
                parse::compress_lazy::<true, false, true>(
                    chains,
                    ctx,
                    seqs(),
                    params,
                    data,
                    &mut w,
                );
            }
            Strategy::NearOptimal(n) => n.compress(ctx, data, &mut w),
        }
        w.finish()
    }
}

/// The most raw DEFLATE data `len` bytes can compress to: stored blocks, at worst, of which
/// there are at most one per [`MIN_BLOCK_LENGTH`] bytes, each with a 5-byte header.
pub(crate) fn compress_bound(len: usize) -> usize {
    5 * len.div_ceil(MIN_BLOCK_LENGTH).max(1) + len
}

/// Writes `data` as stored blocks of up to 65,535 bytes.
fn compress_stored(data: &[u8], out: &mut [u8]) -> Option<usize> {
    let mut w = BitWriter::new(out);
    if data.is_empty() {
        w.write_stored(true, data);
    }
    let mut pieces = data.chunks(usize::from(u16::MAX)).peekable();
    while let Some(piece) = pieces.next() {
        w.write_stored(pieces.peek().is_none(), piece);
        if w.overflowed() {
            return None;
        }
    }
    w.finish()
}
