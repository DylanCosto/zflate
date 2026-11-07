//! Finding earlier copies of the bytes at a position, as libdeflate does: a hash table of
//! the latest positions for level 1, and hash chains for levels 2 to 9.
//!
//! Positions are `i16`, relative to a base that moves up 32 KiB whenever a position reaches
//! it, so the tables stay small; moving it subtracts 32,768 from every entry, and entries
//! that would fall out of the window become [`EMPTY`]. Inputs up to [`SMALL_INPUT_MAX`]
//! bytes, such as BGZF blocks, keep absolute 16-bit positions instead; level 1 then stores
//! each beside a tag that identifies its four bytes.

use crate::alloc::try_boxed_array;

/// DEFLATE's window.
pub(super) const WINDOW_SIZE: usize = 1 << 15;
/// An entry with no position: always out of the window.
const EMPTY: i16 = i16::MIN;

/// Inputs this short have their positions in the tables as they are, with no window to
/// move. The last position must not be [`SMALL_EMPTY`]'s.
pub(super) const SMALL_INPUT_MAX: usize = u16::MAX as usize;
/// An entry with no position in a small input: 65,535, after every position.
const SMALL_EMPTY: i16 = -1;

/// The hash table's and chains' heads' sizes, as libdeflate's.
const TABLE_HASH_ORDER: u32 = 15;
const HASH3_ORDER: u32 = 15;
const HASH4_ORDER: u32 = 16;
/// Where the hash chains' tables start in their allocation. The links come first: following
/// a link is the search's critical path, and at offset 0 it is a mask and a load.
const HC_NEXT: usize = 0;
const HC_HASH3: usize = HC_NEXT + WINDOW_SIZE;
const HC_HASH4: usize = HC_HASH3 + (1 << HASH3_ORDER);
const HC_TABLES_LEN: usize = HC_HASH4 + (1 << HASH4_ORDER);

/// The bytes left a position needs for the hash table to look there: a 4-byte match and the
/// next position's hash.
pub(super) const TABLE_REQUIRED_BYTES: usize = 5;
/// Only look for redundant inserts on long matches, outside the common short-match path.
const MIN_PERIODIC_SKIP: usize = 31;

/// The hashes' multiplier, libdeflate's. It is odd, so multiplying by it is invertible:
/// level 1's short table relies on that to identify four bytes by their product.
const HASH_MULTIPLIER: u32 = 0x1E35_A7BD;
const _: () = assert!(HASH_MULTIPLIER % 2 == 1);

#[inline(always)]
fn hash(seq: u32, bits: u32) -> usize {
    (seq.wrapping_mul(HASH_MULTIPLIER) >> (32 - bits)) as usize
}

/// Asks for the cache line holding `entry`, which the next position will read and write.
#[inline(always)]
fn prefetch<T>(entry: &T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: a prefetch is only a hint; it reads nothing and cannot fault.
    unsafe {
        use std::arch::x86_64::{_MM_HINT_T0, _mm_prefetch};
        _mm_prefetch::<_MM_HINT_T0>(std::ptr::from_ref(entry).cast::<i8>());
    }
    // Not under Miri, which runs no inline assembly; a prefetch changes no result.
    #[cfg(all(target_arch = "aarch64", not(miri)))]
    #[expect(
        clippy::pointers_in_nomem_asm_block,
        reason = "`prfm` takes the address as a hint and reads and writes nothing"
    )]
    // SAFETY: as above: `prfm` is a hint that reads and writes nothing (so `nomem`), and
    // cannot fault.
    // (`core::arch::aarch64::_prefetch` is not stable; libdeflate's C gets this from
    // `__builtin_prefetch`.)
    unsafe {
        std::arch::asm!(
            "prfm pldl1keep, [{0}]",
            in(reg) std::ptr::from_ref(entry),
            options(nomem, nostack, preserves_flags)
        );
    }
    #[cfg(any(miri, not(any(target_arch = "x86_64", target_arch = "aarch64"))))]
    let _ = entry;
}

#[inline(always)]
pub(super) fn load_u32(data: &[u8], pos: usize) -> u32 {
    u32::from_le_bytes(data[pos..pos + 4].try_into().expect("4 bytes"))
}

/// Loads 4 bytes at `p` without a bounds check.
///
/// # Safety
///
/// `p + 4 <= data.len()`.
#[inline(always)]
unsafe fn load_u32_unchecked(data: &[u8], p: usize) -> u32 {
    debug_assert!(p + 4 <= data.len());
    // SAFETY: the caller's.
    // Little-endian, as the 3-byte mask and the hashes assume; free on little-endian CPUs.
    u32::from_le(unsafe { data.as_ptr().add(p).cast::<u32>().read_unaligned() })
}

/// How far the bytes at `pos` and `m` (before it) agree, from `len` up to `max_len`.
#[inline(always)]
pub(super) fn lz_extend(data: &[u8], pos: usize, m: usize, len: u32, max_len: u32) -> u32 {
    let max = max_len as usize;
    let (a, b) = (&data[pos..pos + max], &data[m..m + max]);
    let mut len = len as usize;
    let word = |len: usize| {
        u64::from_le_bytes(a[len..len + 8].try_into().expect("8 bytes"))
            ^ u64::from_le_bytes(b[len..len + 8].try_into().expect("8 bytes"))
    };
    // Matches are often long: with room, compare four words before looping, as libdeflate
    // does.
    if max - len >= 32 {
        for _ in 0..4 {
            let x = word(len);
            if x != 0 {
                return (len + (x.trailing_zeros() >> 3) as usize) as u32;
            }
            len += 8;
        }
    }
    while len + 8 <= max {
        let x = u64::from_le_bytes(a[len..len + 8].try_into().expect("8 bytes"))
            ^ u64::from_le_bytes(b[len..len + 8].try_into().expect("8 bytes"));
        if x != 0 {
            return (len + (x.trailing_zeros() >> 3) as usize) as u32;
        }
        len += 8;
    }
    while len < max && a[len] == b[len] {
        len += 1;
    }
    len as u32
}

/// [`lz_extend`] without bounds checks, for the hash chains' searches.
///
/// # Safety
///
/// `m < pos` and `pos + max_len <= data.len()`.
#[inline(always)]
unsafe fn lz_extend_unchecked(data: &[u8], pos: usize, m: usize, len: u32, max_len: u32) -> u32 {
    debug_assert!(m < pos && pos + max_len as usize <= data.len() && len <= max_len);
    let max = max_len as usize;
    let p = data.as_ptr();
    let mut len = len as usize;
    // SAFETY: each word read is 8 bytes at `len + 8 <= max` from `pos` or from `m`, which
    // is before `pos`, and `pos + max <= data.len()`.
    let word = |len: usize| unsafe {
        // Little-endian, so that the first differing byte is the lowest.
        u64::from_le(
            p.add(pos + len).cast::<u64>().read_unaligned()
                ^ p.add(m + len).cast::<u64>().read_unaligned(),
        )
    };
    if max - len >= 32 {
        for _ in 0..4 {
            let x = word(len);
            if x != 0 {
                return (len + (x.trailing_zeros() >> 3) as usize) as u32;
            }
            len += 8;
        }
    }
    while len + 8 <= max {
        let x = word(len);
        if x != 0 {
            return (len + (x.trailing_zeros() >> 3) as usize) as u32;
        }
        len += 8;
    }
    // SAFETY: as above, a byte at `len < max`.
    while len < max && unsafe { *p.add(pos + len) == *p.add(m + len) } {
        len += 1;
    }
    len as u32
}

/// Moves positions down a window: positive ones by 32,768, the rest out of the window. A
/// saturating add of -32,768 does both, as libdeflate's SSE2 version does: one vector add
/// for 8 entries.
fn rebase(table: &mut [i16]) {
    for entry in table {
        *entry = entry.saturating_add(i16::MIN);
    }
}

/// An empty table, boxed without building it on the stack, or `None` if it cannot be
/// allocated.
fn try_boxed_table<const N: usize>() -> Option<Box<[i16; N]>> {
    try_boxed_array(EMPTY)
}

/// Where the window's position 0 is in the input.
struct Window {
    base: usize,
}

impl Window {
    /// The position of `pos` in the window, moving the window (and the tables, with
    /// `rebase`) up when `pos` reaches its end.
    #[inline(always)]
    fn pos(&mut self, pos: usize, rebase: impl FnOnce()) -> usize {
        let mut cur = pos - self.base;
        if cur == WINDOW_SIZE {
            rebase();
            self.base += WINDOW_SIZE;
            cur = 0;
        }
        cur
    }

    /// The input position of the window position `node`, which may be in the window before.
    #[inline(always)]
    fn at(&self, node: i16) -> usize {
        self.base.wrapping_add_signed(isize::from(node))
    }

    /// [`at`](Self::at), if it is before `pos`. Positions come from the tables, which only
    /// hold earlier ones; checking costs a branch that is never taken, and makes the
    /// unchecked loads at the position sound whatever the tables hold.
    #[inline(always)]
    fn before(&self, node: i16, pos: usize) -> Option<usize> {
        let m = self.at(node);
        (m < pos).then_some(m)
    }
}

/// Positions at or below this are out of the window from `cur`.
#[inline(always)]
fn cutoff(cur: usize) -> i16 {
    (cur as i32 - WINDOW_SIZE as i32) as i16
}

/// Level 1's table for inputs up to [`SMALL_INPUT_MAX`] bytes: 2^16 buckets of two entries,
/// each a 16-bit tag above a 16-bit absolute position (512 KiB). An empty entry, `u32::MAX`,
/// carries tag 0xFFFF and so passes the tag test for some keys, as position 65,535; it is
/// rejected because that position is never before the current one (no position of a short
/// input reaches it), so its offset is out of the window. Keep the offset test.
type ShortBuckets = [[u32; 2]; 1 << 16];
/// Level 1's table for longer inputs, with `i16` positions in the sliding window, at the
/// start of the short table's allocation (128 KiB).
type RollingBuckets = [[i16; 2]; 1 << TABLE_HASH_ORDER];

// The rolling table fits in the short table's allocation, at its alignment.
const _: () = assert!(
    size_of::<RollingBuckets>() <= size_of::<ShortBuckets>()
        && align_of::<RollingBuckets>() <= align_of::<ShortBuckets>()
);

/// Level 1's two-candidate match finder. Short inputs use the upper 16 bits of
/// the hash to select a bucket and store the lower 16 bits alongside each absolute
/// position. The odd multiplier is invertible, so the index and tag together
/// identify the four bytes exactly.
///
/// The short table is four times libdeflate's, and clearing it for each input is about 3%
/// of level 1's time on a BAM; it no longer fits in a 128 KiB L1 cache, as libdeflate's
/// does on Apple silicon.
pub(super) struct HashTable {
    buckets: Box<ShortBuckets>,
    window: Window,
    next_hash: usize,
}

impl HashTable {
    pub(super) fn try_new() -> Option<Self> {
        Some(Self {
            buckets: try_boxed_array([u32::MAX; 2])?,
            window: Window { base: 0 },
            next_hash: 0,
        })
    }

    /// The rolling table uses the start of the short table's allocation.
    fn rolling_buckets(buckets: &mut ShortBuckets) -> &mut RollingBuckets {
        // SAFETY: the allocation is initialized, and at least as large and as aligned as
        // the rolling table (asserted above). Every bit pattern is valid, and the
        // returned borrow excludes access in the other representation.
        unsafe { &mut *buckets.as_mut_ptr().cast::<RollingBuckets>() }
    }

    /// Empties the table and returns the short input's first key. The parser
    /// carries this key in a local to keep it in a register.
    pub(super) fn reset<const SMALL: bool>(&mut self, data: &[u8]) -> u32 {
        self.window.base = 0;
        let first = data
            .first_chunk::<4>()
            .map_or(0, |word| u32::from_le_bytes(*word));
        if SMALL {
            self.buckets.fill([u32::MAX; 2]);
            Self::key(first)
        } else {
            Self::rolling_buckets(&mut self.buckets).fill([EMPTY; 2]);
            self.next_hash = hash(first, TABLE_HASH_ORDER);
            0
        }
    }

    #[inline(always)]
    fn key(word: u32) -> u32 {
        word.wrapping_mul(HASH_MULTIPLIER)
    }

    /// # Safety
    ///
    /// There are at least max_len bytes at pos, with 5 <= max_len <= 258.
    /// The caller has reset the short table for this input and visits positions
    /// in order. The input fits in 16-bit positions and nice_len <= max_len.
    #[inline(always)]
    unsafe fn longest_match_small(
        &mut self,
        next_key: &mut u32,
        data: &[u8],
        pos: usize,
        max_len: u32,
        nice_len: u32,
    ) -> (u32, u32) {
        let key = *next_key;
        // SAFETY: there are at least five bytes at pos.
        *next_key = Self::key(unsafe { load_u32_unchecked(data, pos + 1) });
        let slot = &mut self.buckets[(key >> 16) as usize];
        let first = slot[0];
        let second = slot[1];
        let tag = key << 16;
        // Store the two entries independently. Consecutive positions with the
        // same key can forward the new entry without waiting on the old one.
        slot[0] = tag | pos as u32;
        slot[1] = first;
        let hit = (first ^ tag) >> 16 == 0;
        let first = if hit { first } else { second };
        let offset = pos.wrapping_sub(usize::from(first as u16));
        if (first ^ tag) >> 16 != 0 || offset.wrapping_sub(1) >= WINDOW_SIZE - 1 {
            return (0, 0);
        }
        // SAFETY: the full key identifies four equal bytes. The offset is
        // positive and at most pos, and max_len bytes are available at pos.
        let best_len = unsafe { lz_extend_unchecked(data, pos, pos - offset, 4, max_len) };
        if !hit || best_len >= nice_len || (second ^ tag) >> 16 != 0 {
            return (best_len, offset as u32);
        }
        let m = usize::from(second as u16);
        let second_offset = pos.wrapping_sub(m);
        if second_offset.wrapping_sub(1) < WINDOW_SIZE - 1 {
            let probe = best_len as usize - 3;
            // SAFETY: m precedes pos, and best_len is below nice_len <= max_len.
            if unsafe {
                load_u32_unchecked(data, m + probe) == load_u32_unchecked(data, pos + probe)
            } {
                // SAFETY: the tag identifies four equal bytes. The candidate
                // precedes pos and max_len bytes fit at pos.
                let len = unsafe { lz_extend_unchecked(data, pos, m, 4, max_len) };
                if len > best_len {
                    return (len, second_offset as u32);
                }
            }
        }
        (best_len, offset as u32)
    }

    fn rebase(table: &mut [[i16; 2]]) {
        rebase(table.as_flattened_mut());
    }

    /// Inserts `pos` and returns the longest match there (length and offset), or length 0.
    ///
    /// # Safety
    ///
    /// There are at least max_len bytes at pos, with 5 <= max_len <= 258 and
    /// nice_len <= max_len. The table was reset for this input in the same mode.
    /// In SMALL mode the input fits in 16-bit positions.
    #[inline(always)]
    pub(super) unsafe fn longest_match<const SMALL: bool>(
        &mut self,
        next_key: &mut u32,
        data: &[u8],
        pos: usize,
        max_len: u32,
        nice_len: u32,
    ) -> (u32, u32) {
        if SMALL {
            // SAFETY: the caller guarantees the short matcher's requirements.
            return unsafe { self.longest_match_small(next_key, data, pos, max_len, nice_len) };
        }
        let table = Self::rolling_buckets(&mut self.buckets);
        let cur = self.window.pos(pos, || Self::rebase(table));
        let cutoff = cutoff(cur);
        let h = self.next_hash & ((1 << TABLE_HASH_ORDER) - 1);
        self.next_hash = hash(load_u32(data, pos + 1), TABLE_HASH_ORDER);
        prefetch(&table[self.next_hash & ((1 << TABLE_HASH_ORDER) - 1)]);
        let seq = load_u32(data, pos);
        let bucket = &mut table[h];
        let first = bucket[0];
        bucket[0] = cur as i16;
        if first <= cutoff {
            return (0, 0);
        }
        let second = bucket[1];
        bucket[1] = first;
        let mut best_len = 0;
        let mut best = pos;
        let m = self.window.at(first);
        if load_u32(data, m) == seq {
            best_len = lz_extend(data, pos, m, 4, max_len);
            best = m;
            if second <= cutoff || best_len >= nice_len {
                return (best_len, (pos - best) as u32);
            }
            let m = self.window.at(second);
            if load_u32(data, m) == seq
                && load_u32(data, m + best_len as usize - 3)
                    == load_u32(data, pos + best_len as usize - 3)
            {
                let len = lz_extend(data, pos, m, 4, max_len);
                if len > best_len {
                    best_len = len;
                    best = m;
                }
            }
        } else {
            if second <= cutoff {
                return (0, 0);
            }
            let m = self.window.at(second);
            if load_u32(data, m) == seq {
                best_len = lz_extend(data, pos, m, 4, max_len);
                best = m;
            }
        }
        (best_len, (pos - best) as u32)
    }

    /// Inserts the `count` positions from `pos` (the one after a match's start), without
    /// looking for matches; does nothing near the end, where no more matches are sought.
    /// The copied bytes repeat every `period` bytes, the match's offset.
    #[inline(always)]
    pub(super) fn skip<const SMALL: bool>(
        &mut self,
        next_key: &mut u32,
        data: &[u8],
        pos: usize,
        count: usize,
        period: u32,
    ) {
        if data.len().saturating_sub(pos).saturating_sub(count) < TABLE_REQUIRED_BYTES {
            crate::cold_path();
            return;
        }
        if SMALL {
            let end = pos + count;
            // A copied run repeats every period bytes. The final two periods
            // retain both candidates for every key before the three trailing
            // positions whose keys can extend beyond the run.
            let keep = 2 * period as usize + 3;
            let start = if count >= MIN_PERIODIC_SKIP {
                end - count.min(keep)
            } else {
                pos
            };
            let mut key = if start == pos {
                *next_key
            } else {
                // SAFETY: start <= end, and the guard leaves five bytes at end.
                Self::key(unsafe { load_u32_unchecked(data, start) })
            };
            for p in start..end {
                let slot = &mut self.buckets[(key >> 16) as usize];
                slot[1] = slot[0];
                slot[0] = (key << 16) | p as u32;
                // SAFETY: the skip guard leaves five bytes at pos + count.
                key = Self::key(unsafe { load_u32_unchecked(data, p + 1) });
            }
            *next_key = key;
            return;
        }
        let table = Self::rolling_buckets(&mut self.buckets);
        // Move the window first if the positions reach its end; the first of them may then
        // be negative, in the window before.
        let mut cur = (pos - self.window.base) as isize;
        if cur as usize + count > WINDOW_SIZE {
            Self::rebase(table);
            self.window.base += WINDOW_SIZE;
            cur -= WINDOW_SIZE as isize;
        }
        let mut h = self.next_hash;
        for p in pos..pos + count {
            let bucket = &mut table[h];
            bucket[1] = bucket[0];
            bucket[0] = cur as i16;
            cur += 1;
            h = hash(load_u32(data, p + 1), TABLE_HASH_ORDER);
        }
        prefetch(&table[h]);
        self.next_hash = h;
    }
}

/// Levels 2 to 9's match finder: the latest position of each 3-byte sequence's hash, and
/// chains through all positions of each 4-byte sequence's hash.
///
/// Its methods take `FIND3`: whether the parser looks for 3-byte matches. Without them,
/// the 3-byte heads are neither kept nor read.
pub(super) struct HashChains {
    /// Each position's predecessor in its 4-byte chain, the 3-byte heads, and the 4-byte
    /// heads, in one allocation, so the search needs one base pointer and not three: it
    /// is short of registers.
    tables: Box<[i16; HC_TABLES_LEN]>,
    window: Window,
}

/// What a match search found: the longest match's length and offset. With no match longer
/// than the search started with, the length is that and the offset 0.
#[derive(Clone, Copy, Debug)]
pub(super) struct Match {
    pub(super) len: u32,
    pub(super) offset: u32,
}

impl HashChains {
    pub(super) fn try_new() -> Option<Self> {
        Some(Self {
            tables: try_boxed_table()?,
            window: Window { base: 0 },
        })
    }

    /// Empties the tables for a new input, `SMALL` if it is at most [`SMALL_INPUT_MAX`]
    /// bytes. The chains need no clearing: a position's link is written when it is
    /// inserted, before anything can follow it. Without `FIND3`, neither do the 3-byte
    /// heads, which are not used.
    /// Returns the first position's 3- and 4-byte hashes. The parser keeps the next
    /// position's hashes, and passes them to each search and skip: a local of the parser
    /// stays in a register, where a field here would be stored and reloaded each search.
    pub(super) fn reset<const SMALL: bool, const FIND3: bool>(
        &mut self,
        data: &[u8],
    ) -> [usize; 2] {
        debug_assert_eq!(SMALL, data.len() <= SMALL_INPUT_MAX);
        let heads = if FIND3 { HC_HASH3 } else { HC_HASH4 };
        self.tables[heads..].fill(if SMALL { SMALL_EMPTY } else { EMPTY });
        self.window.base = 0;
        if data.len() >= 4 {
            Self::hashes::<FIND3>(load_u32(data, 0))
        } else {
            [0; 2]
        }
    }

    /// The 3-byte (with `FIND3`; else 0) and 4-byte hashes of `seq`'s bytes.
    #[inline(always)]
    fn hashes<const FIND3: bool>(seq: u32) -> [usize; 2] {
        let h3 = if FIND3 {
            hash(seq & 0xff_ffff, HASH3_ORDER)
        } else {
            0
        };
        [h3, hash(seq, HASH4_ORDER)]
    }

    #[inline(always)]
    fn hash3_slot(h: usize) -> usize {
        HC_HASH3 + (h & ((1 << HASH3_ORDER) - 1))
    }

    #[inline(always)]
    fn hash4_slot(h: usize) -> usize {
        HC_HASH4 + (h & ((1 << HASH4_ORDER) - 1))
    }

    /// The chain link of the window position `cur`, or of `node`, which may be in the
    /// window before.
    #[inline(always)]
    fn next_slot(cur: usize) -> usize {
        HC_NEXT + (cur & (WINDOW_SIZE - 1))
    }

    /// The position `pos` as the tables store it: for a small input the position itself,
    /// which needs no window; otherwise its place in the window, moving the window (and
    /// the tables) when `pos` reaches its end.
    #[inline(always)]
    fn cur<const SMALL: bool, const FIND3: bool>(&mut self, pos: usize) -> usize {
        if SMALL {
            pos
        } else {
            let tables = &mut self.tables;
            self.window.pos(pos, || {
                // Without `FIND3`, the 3-byte heads are not used.
                let (links_and_heads3, heads4) = tables.split_at_mut(HC_HASH4);
                rebase(if FIND3 {
                    links_and_heads3
                } else {
                    &mut links_and_heads3[..HC_HASH3]
                });
                rebase(heads4);
            })
        }
    }

    /// The input position of the table entry `node`, if it is a candidate for a match at
    /// `pos`: in the window, and before `pos`. For a small input, one unsigned compare
    /// checks both, and rejects [`SMALL_EMPTY`] too. Entries are always earlier positions,
    /// but checking costs little and makes the unchecked loads at the position sound
    /// whatever the tables hold.
    #[inline(always)]
    fn candidate<const SMALL: bool>(&self, node: i16, pos: usize, cutoff: i16) -> Option<usize> {
        if SMALL {
            let m = usize::from(node as u16);
            (pos.wrapping_sub(m).wrapping_sub(1) < WINDOW_SIZE - 1).then_some(m)
        } else if node > cutoff {
            self.window.before(node, pos)
        } else {
            None
        }
    }

    /// Inserts `pos` and returns the longest match there longer than `best_len`, looking at
    /// most `max_depth` positions down the chain and stopping at `nice_len`. Fewer than 5
    /// bytes from the end, it neither inserts nor looks. Without `FIND3`, `best_len` is at
    /// least 3.
    #[expect(
        clippy::too_many_arguments,
        reason = "a search's hash state, input, position, threshold and limits, all inlined"
    )]
    #[inline(always)]
    pub(super) fn longest_match<const SMALL: bool, const FIND3: bool>(
        &mut self,
        next_hashes: &mut [usize; 2],
        data: &[u8],
        pos: usize,
        mut best_len: u32,
        max_len: u32,
        nice_len: u32,
        max_depth: u32,
    ) -> Match {
        let cur = self.cur::<SMALL, FIND3>(pos);
        let cutoff = cutoff(cur);
        let none = Match {
            len: best_len,
            offset: 0,
        };
        if max_len < 5 {
            crate::cold_path();
            return none;
        }
        // The candidates' loads below read at most `max_len` bytes from `pos`. The parsers
        // never ask for more than there is; returning, not panicking, keeps the search free
        // of a panic path.
        if pos + max_len as usize > data.len() {
            crate::cold_path();
            return none;
        }
        // The loads at the end of the best match need `best_len < max_len`, which the search
        // keeps by stopping at `nice_len`: kept here, it holds whatever the parser passes.
        let nice_len = nice_len.min(max_len);
        debug_assert!(FIND3 || best_len >= 3);
        let [h3, h4] = *next_hashes;
        let t = &mut self.tables;
        let node3 = if FIND3 {
            t[Self::hash3_slot(h3)]
        } else {
            EMPTY
        };
        let mut node4 = t[Self::hash4_slot(h4)];
        if FIND3 {
            t[Self::hash3_slot(h3)] = cur as i16;
        }
        t[Self::hash4_slot(h4)] = cur as i16;
        t[Self::next_slot(cur)] = node4;
        // SAFETY (for every unchecked load and extension below): the check above puts
        // `max_len` bytes, at least 5, from `pos` in `data`; candidates are before `pos`
        // (`candidate`), and a load at `end` is within `max_len`, as `best_len < nice_len
        // <= max_len` (clamped above).
        let [n3, n4] = Self::hashes::<FIND3>(unsafe { load_u32_unchecked(data, pos + 1) });
        *next_hashes = [n3, n4];
        if FIND3 {
            prefetch(&t[Self::hash3_slot(n3)]);
        }
        prefetch(&t[Self::hash4_slot(n4)]);

        let mut depth = max_depth;
        let mut best = pos;
        let found = |best_len: u32, best: usize| Match {
            len: best_len,
            offset: (pos - best) as u32,
        };
        let link = |node: i16| self.tables[Self::next_slot(usize::from(node as u16))];
        if best_len < 4 {
            let seq4 = unsafe { load_u32_unchecked(data, pos) };
            if FIND3 {
                // With no position of the same 3-byte hash in the window, there is no match.
                let Some(m3) = self.candidate::<SMALL>(node3, pos, cutoff) else {
                    return none;
                };
                if best_len < 3
                    && unsafe { load_u32_unchecked(data, m3) } & 0xff_ffff == seq4 & 0xff_ffff
                {
                    best_len = 3;
                    best = m3;
                }
            }
            // Find the first 4-byte match.
            let m = loop {
                let Some(m) = self.candidate::<SMALL>(node4, pos, cutoff) else {
                    return found(best_len, best);
                };
                if unsafe { load_u32_unchecked(data, m) } == seq4 {
                    break m;
                }
                node4 = link(node4);
                depth -= 1;
                if depth == 0 {
                    return found(best_len, best);
                }
            };
            best = m;
            best_len = unsafe { lz_extend_unchecked(data, pos, m, 4, max_len) };
            if best_len >= nice_len {
                return found(best_len, best);
            }
            node4 = link(node4);
            depth -= 1;
            if depth == 0 {
                return found(best_len, best);
            }
        } else if best_len >= nice_len {
            return none;
        }
        loop {
            // Skip to a candidate that could be longer: it must agree at the 4 bytes ending
            // one past the best length, and at its first 4.
            let end = best_len as usize - 3;
            let (seq_end, seq4) = unsafe {
                (
                    load_u32_unchecked(data, pos + end),
                    load_u32_unchecked(data, pos),
                )
            };
            let m = loop {
                let Some(m) = self.candidate::<SMALL>(node4, pos, cutoff) else {
                    return found(best_len, best);
                };
                if unsafe { load_u32_unchecked(data, m + end) } == seq_end
                    && unsafe { load_u32_unchecked(data, m) } == seq4
                {
                    break m;
                }
                node4 = link(node4);
                depth -= 1;
                if depth == 0 {
                    return found(best_len, best);
                }
            };
            let len = unsafe { lz_extend_unchecked(data, pos, m, 4, max_len) };
            if len > best_len {
                best_len = len;
                best = m;
                if best_len >= nice_len {
                    break;
                }
            }
            node4 = link(node4);
            depth -= 1;
            if depth == 0 {
                break;
            }
        }
        found(best_len, best)
    }

    /// Inserts the remaining positions of a match. At level 5, no search
    /// can reach more than depth entries in a chain. In a periodic run
    /// every key in the omitted prefix has at least depth newer entries.
    /// Keep three more positions whose four-byte keys can extend past the run.
    /// `depth` must be the deepest any search of this input goes, or the output
    /// changes (still valid, but with matches missed).
    #[inline(always)]
    pub(super) fn skip_match<const SMALL: bool, const FIND3: bool>(
        &mut self,
        next_hashes: &mut [usize; 2],
        data: &[u8],
        mut pos: usize,
        mut count: usize,
        period: u32,
        depth: u32,
    ) {
        if SMALL && !FIND3 && count >= MIN_PERIODIC_SKIP {
            let keep = depth as usize * period as usize + 3;
            if count > keep && data.len() - pos >= count + 5 {
                pos += count - keep;
                count = keep;
                *next_hashes = Self::hashes::<FIND3>(load_u32(data, pos));
            }
        }
        self.skip::<SMALL, FIND3>(next_hashes, data, pos, count);
    }

    /// Inserts the `count` positions from `pos` without looking for matches; does nothing
    /// near the end, where no more matches are sought.
    #[inline(always)]
    pub(super) fn skip<const SMALL: bool, const FIND3: bool>(
        &mut self,
        next_hashes: &mut [usize; 2],
        data: &[u8],
        pos: usize,
        count: usize,
    ) {
        if count + 5 > data.len() - pos {
            crate::cold_path();
            return;
        }
        let [mut h3, mut h4] = *next_hashes;
        for p in pos..pos + count {
            let cur = self.cur::<SMALL, FIND3>(p);
            let t = &mut self.tables;
            if FIND3 {
                t[Self::hash3_slot(h3)] = cur as i16;
            }
            t[Self::next_slot(cur)] = t[Self::hash4_slot(h4)];
            t[Self::hash4_slot(h4)] = cur as i16;
            // SAFETY: `count + 5 <= data.len() - pos`, checked above, so the 4 bytes at
            // `p + 1 < pos + count` are in `data`.
            [h3, h4] = Self::hashes::<FIND3>(unsafe { load_u32_unchecked(data, p + 1) });
        }
        if FIND3 {
            prefetch(&self.tables[Self::hash3_slot(h3)]);
        }
        prefetch(&self.tables[Self::hash4_slot(h4)]);
        *next_hashes = [h3, h4];
    }
}

/// The binary trees' hash orders, as libdeflate's: two ways for 3-byte sequences.
const BT_HASH3_ORDER: u32 = 16;
const BT_HASH4_ORDER: u32 = 16;

/// The bytes left a position needs for the binary trees to look there.
pub(super) const TREES_REQUIRED_BYTES: u32 = 5;

/// A match the near-optimal parser may use: its length and offset. In the match cache, each
/// position's matches are followed by an entry with their count as `length` and the
/// position's byte as `offset`.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct LzMatch {
    pub(super) length: u16,
    pub(super) offset: u16,
}

/// Levels 10 to 12's match finder: for each 4-byte sequence's hash, a binary tree of the
/// positions with it, ordered by the bytes that follow, which yields every match that is
/// longer than all the matches found before it; and the latest two positions of each 3-byte
/// sequence's hash, for 3-byte matches.
pub(super) struct BinaryTrees {
    hash3: Box<[[i16; 2]; 1 << BT_HASH3_ORDER]>,
    hash4: Box<[i16; 1 << BT_HASH4_ORDER]>,
    /// Each position's two children: the positions whose following bytes compare lower and
    /// higher.
    children: Box<[i16; 2 * WINDOW_SIZE]>,
    window: Window,
    next_hashes: [usize; 2],
}

impl BinaryTrees {
    pub(super) fn try_new() -> Option<Self> {
        Some(Self {
            hash3: try_boxed_array([EMPTY; 2])?,
            hash4: try_boxed_table()?,
            children: try_boxed_table()?,
            window: Window { base: 0 },
            next_hashes: [0; 2],
        })
    }

    /// Empties the tables for a new input. The children need no clearing: a position's are
    /// written when it is inserted.
    pub(super) fn reset(&mut self, data: &[u8]) {
        self.hash3.fill([EMPTY; 2]);
        self.hash4.fill(EMPTY);
        self.window.base = 0;
        self.next_hashes = if data.len() >= 4 {
            Self::hashes(load_u32(data, 0))
        } else {
            [0; 2]
        };
    }

    #[inline(always)]
    fn hashes(seq: u32) -> [usize; 2] {
        [
            hash(seq & 0xff_ffff, BT_HASH3_ORDER),
            hash(seq, BT_HASH4_ORDER),
        ]
    }

    /// Inserts `pos` and appends to `matches` each match there longer than those before it,
    /// up to `max_len`, looking at most `max_depth` nodes down the tree and stopping at
    /// `nice_len`. `max_len` must be at least [`TREES_REQUIRED_BYTES`].
    #[inline(always)]
    pub(super) fn get_matches(
        &mut self,
        data: &[u8],
        pos: usize,
        max_len: u32,
        nice_len: u32,
        max_depth: u32,
        matches: &mut Vec<LzMatch>,
    ) {
        self.advance::<true>(data, pos, max_len, nice_len, max_depth, matches);
    }

    /// Inserts `pos` without recording matches.
    #[inline(always)]
    pub(super) fn skip(&mut self, data: &[u8], pos: usize, nice_len: u32, max_depth: u32) {
        // Nothing is recorded, so an empty vector, which does not allocate, will do.
        self.advance::<false>(data, pos, nice_len, nice_len, max_depth, &mut Vec::new());
    }

    /// Inserts `pos` at the root of its tree: walks down from the old root, comparing the
    /// bytes at each node with those at `pos`, and splits the tree into the nodes lower
    /// and higher than `pos`, which become its children. The walk passes the longest
    /// matches on the way, which `RECORD` records.
    #[inline(always)]
    fn advance<const RECORD: bool>(
        &mut self,
        data: &[u8],
        pos: usize,
        max_len: u32,
        nice_len: u32,
        max_depth: u32,
        matches: &mut Vec<LzMatch>,
    ) {
        debug_assert!(max_len >= TREES_REQUIRED_BYTES && pos + max_len as usize <= data.len());
        let (hash3, hash4, children) = (&mut self.hash3, &mut self.hash4, &mut self.children);
        let cur = self.window.pos(pos, || {
            rebase(hash3.as_flattened_mut());
            rebase(&mut hash4[..]);
            rebase(&mut children[..]);
        });
        let cutoff = cutoff(cur);
        let [h3, h4] = self.next_hashes;
        let (h3, h4) = (
            h3 & ((1 << BT_HASH3_ORDER) - 1),
            h4 & ((1 << BT_HASH4_ORDER) - 1),
        );
        self.next_hashes = Self::hashes(load_u32(data, pos + 1));
        prefetch(&self.hash3[self.next_hashes[0] & ((1 << BT_HASH3_ORDER) - 1)]);
        prefetch(&self.hash4[self.next_hashes[1] & ((1 << BT_HASH4_ORDER) - 1)]);

        let [node3, node3_2] = self.hash3[h3];
        self.hash3[h3] = [cur as i16, node3];
        if RECORD && node3 > cutoff {
            let seq3 = load_u32(data, pos) & 0xff_ffff;
            let m = self.window.at(node3);
            if load_u32(data, m) & 0xff_ffff == seq3 {
                matches.push(LzMatch {
                    length: 3,
                    offset: (pos - m) as u16,
                });
            } else if node3_2 > cutoff {
                let m = self.window.at(node3_2);
                if load_u32(data, m) & 0xff_ffff == seq3 {
                    matches.push(LzMatch {
                        length: 3,
                        offset: (pos - m) as u16,
                    });
                }
            }
        }

        let mut node = self.hash4[h4];
        self.hash4[h4] = cur as i16;
        // Where the next lower and higher nodes go: at first, `pos`'s own children.
        let (mut pending_lt, mut pending_gt) = (2 * cur, 2 * cur + 1);
        if node <= cutoff {
            self.children[pending_lt] = EMPTY;
            self.children[pending_gt] = EMPTY;
            return;
        }
        let left = |node: i16| 2 * (usize::from(node as u16) & (WINDOW_SIZE - 1));
        let (mut best_lt_len, mut best_gt_len) = (0, 0);
        let mut len = 0;
        let mut best_len = 3;
        let mut depth = max_depth;
        loop {
            let m = self.window.at(node);
            // The node agrees with `pos` for at least `len` bytes: the shorter of the
            // lower and higher bounds' agreements.
            if data[m + len as usize] == data[pos + len as usize] {
                len = lz_extend(data, pos, m, len + 1, max_len);
                if !RECORD || len > best_len {
                    if RECORD {
                        best_len = len;
                        matches.push(LzMatch {
                            length: len as u16,
                            offset: (pos - m) as u16,
                        });
                    }
                    if len >= nice_len {
                        // Long enough: `pos` takes over the node's children.
                        self.children[pending_lt] = self.children[left(node)];
                        self.children[pending_gt] = self.children[left(node) + 1];
                        return;
                    }
                }
            }
            if data[m + len as usize] < data[pos + len as usize] {
                self.children[pending_lt] = node;
                pending_lt = left(node) + 1;
                node = self.children[pending_lt];
                best_lt_len = len;
                len = len.min(best_gt_len);
            } else {
                self.children[pending_gt] = node;
                pending_gt = left(node);
                node = self.children[pending_gt];
                best_gt_len = len;
                len = len.min(best_lt_len);
            }
            depth -= 1;
            if node <= cutoff || depth == 0 {
                self.children[pending_lt] = EMPTY;
                self.children[pending_gt] = EMPTY;
                return;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Level 5's search depth, the deepest any of its searches goes.
    const DEPTH: u32 = 12;

    /// A short input: random bytes over a small alphabet, then `run` bytes that repeat the
    /// `period` bytes before them, then random bytes, the same pattern again, and more random
    /// bytes, so that later searches find positions inside the run. Returns the input and
    /// where the repeat starts.
    fn periodic(period: usize, run: usize, seed: u64) -> (Vec<u8>, usize) {
        let mut x = seed;
        let mut byte = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[(x % 4) as usize]
        };
        let mut data: Vec<u8> = (0..2000).map(|_| byte()).collect();
        let start = data.len();
        for i in start..start + run {
            data.push(data[i - period]);
        }
        let again = data[start..start + run.min(300)].to_vec();
        data.extend((0..200).map(|_| byte()));
        data.extend(again);
        data.extend((0..100).map(|_| byte()));
        (data, start)
    }

    /// Inserts `count` positions from `pos` into level 1's short table, in pieces too short
    /// to leave any out.
    fn insert_all(table: &mut HashTable, key: &mut u32, data: &[u8], pos: usize, count: usize) {
        for at in (pos..pos + count).step_by(MIN_PERIODIC_SKIP - 1) {
            let n = (MIN_PERIODIC_SKIP - 1).min(pos + count - at);
            table.skip::<true>(key, data, at, n, 1);
        }
    }

    /// Leaving out the inserts of a periodic match ends with level 1's table as inserting
    /// every position would.
    #[test]
    fn level_1_skips_only_inserts_that_cannot_matter() {
        let mut seed = 0x9E37_79B9_7F4A_7C15;
        for period in [1, 2, 3, 4, 7, 13, 64, 129, 300] {
            for len in [31, 32, 33, 40, 100, 200, 258] {
                seed += 1;
                let (data, start) = periodic(period, 258, seed);
                let mut all = HashTable::try_new().expect("table");
                let mut some = HashTable::try_new().expect("table");
                let mut all_key = all.reset::<true>(&data);
                let mut some_key = some.reset::<true>(&data);
                // The positions before the match and its first, as the parser inserts them.
                insert_all(&mut all, &mut all_key, &data, 0, start + 1);
                insert_all(&mut some, &mut some_key, &data, 0, start + 1);
                insert_all(&mut all, &mut all_key, &data, start + 1, len - 1);
                some.skip::<true>(&mut some_key, &data, start + 1, len - 1, period as u32);
                assert!(
                    *all.buckets == *some.buckets && all_key == some_key,
                    "period {period}, length {len}"
                );
            }
        }
    }

    /// Leaving out the inserts of a periodic match at level 5 changes no later search.
    #[test]
    fn level_5_skips_only_inserts_that_cannot_matter() {
        let mut seed = 0x2545_F491_4F6C_DD1D;
        for period in [1, 2, 3, 4, 7, 13, 64, 129, 300] {
            for len in [31, 32, 40, 100, 200, 258] {
                seed += 1;
                let (data, start) = periodic(period, 258, seed);
                let mut all = HashChains::try_new().expect("chains");
                let mut some = HashChains::try_new().expect("chains");
                let mut all_hashes = all.reset::<true, false>(&data);
                let mut some_hashes = some.reset::<true, false>(&data);
                all.skip::<true, false>(&mut all_hashes, &data, 0, start + 1);
                some.skip::<true, false>(&mut some_hashes, &data, 0, start + 1);
                all.skip::<true, false>(&mut all_hashes, &data, start + 1, len - 1);
                some.skip_match::<true, false>(
                    &mut some_hashes,
                    &data,
                    start + 1,
                    len - 1,
                    period as u32,
                    DEPTH,
                );
                // Every later position, searched as deep as level 5's search and as its
                // look-ahead (3), in turn.
                for pos in start + len..data.len() - 5 {
                    let max_len = (data.len() - pos).min(258) as u32;
                    let depth = if pos % 2 == 0 { DEPTH } else { 3 };
                    let a = all.longest_match::<true, false>(
                        &mut all_hashes,
                        &data,
                        pos,
                        3,
                        max_len,
                        max_len,
                        depth,
                    );
                    let s = some.longest_match::<true, false>(
                        &mut some_hashes,
                        &data,
                        pos,
                        3,
                        max_len,
                        max_len,
                        depth,
                    );
                    assert_eq!(
                        (a.len, a.offset),
                        (s.len, s.offset),
                        "period {period}, length {len}, position {pos}"
                    );
                }
            }
        }
    }
}
