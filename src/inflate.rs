//! Raw DEFLATE decompression (RFC 1951) into a caller's buffer, after libdeflate's
//! `decompress_template.h` (MIT).
//!
//! This is the plain version: one symbol per refill, matches copied a byte at a time, every
//! access checked. It accepts and rejects exactly what libdeflate does, so the fast paths can
//! be added one at a time against it.

mod batch;
mod bits;
mod table;

pub(crate) use self::batch::inflate_batch;
use self::bits::{Bits, CONSUMABLE, mask};
use self::table::{
    END_OF_BLOCK, EXCEPTIONAL, LITERAL, LITLEN_ENOUGH, LITLEN_RESULTS, LITLEN_TABLEBITS,
    LITLEN_TABLEBITS_PLAIN, MAX_CODEWORD_LEN, MAX_PRE_CODEWORD_LEN, NUM_LITLEN_SYMS,
    NUM_OFFSET_SYMS, NUM_PRECODE_SYMS, OFFSET_ENOUGH, OFFSET_RESULTS, OFFSET_TABLEBITS, PAIR,
    PRECODE_ENOUGH, PRECODE_RESULTS, PRECODE_TABLEBITS, Pairs, SUBTABLE_POINTER, TableBits,
};
use crate::{DecompressionError, alloc::try_box};

/// The most code lengths a run can write past the last one: a run of 138 zeros with one
/// length left. Sizing `lens` for it saves checking each run before writing it.
const MAX_LENS_OVERRUN: usize = 137;

/// The order of the precode lengths in a dynamic block header (RFC 1951 3.2.7).
const PRECODE_ORDER: [usize; NUM_PRECODE_SYMS] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];

/// The decode tables and their scratch space, kept between calls: too big for the stack,
/// and the static codes' tables can be reused.
#[derive(Debug)]
pub(crate) struct Tables {
    precode: [u32; PRECODE_ENOUGH],
    litlen: [u32; LITLEN_ENOUGH],
    offset: [u32; OFFSET_ENOUGH],
    precode_lens: [u8; NUM_PRECODE_SYMS],
    lens: [u8; NUM_LITLEN_SYMS + NUM_OFFSET_SYMS + MAX_LENS_OVERRUN],
    sorted: [u16; NUM_LITLEN_SYMS],
    litlen_bits: u32,
    /// `litlen` and `offset` hold the static codes' tables.
    static_loaded: bool,
    /// Whether to add [`PAIR`] entries to the literal/length tables; see [`pairs_pay`].
    pairs: bool,
}

impl Tables {
    /// Empty tables, or `None` if they cannot be allocated.
    pub(crate) fn try_new() -> Option<Box<Self>> {
        try_box(Self {
            precode: [0; PRECODE_ENOUGH],
            litlen: [0; LITLEN_ENOUGH],
            offset: [0; OFFSET_ENOUGH],
            precode_lens: [0; NUM_PRECODE_SYMS],
            lens: [0; NUM_LITLEN_SYMS + NUM_OFFSET_SYMS + MAX_LENS_OVERRUN],
            sorted: [0; NUM_LITLEN_SYMS],
            litlen_bits: 0,
            static_loaded: false,
            pairs: false,
        })
    }
}

/// What an inflation produced.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct Inflated {
    /// The input bytes the stream took, up to and including its last block.
    pub(crate) consumed: usize,
    /// The bytes written to the start of the output.
    pub(crate) written: usize,
}

/// Inflates the DEFLATE stream at the start of `data` into `out`.
pub(crate) fn inflate(
    t: &mut Tables,
    data: &[u8],
    out: &mut [u8],
) -> Result<Inflated, DecompressionError> {
    // Decided here, once per call, rather than in the decoder: code added to the decoder's
    // hot function costs more than this saves.
    let pairs = pairs_pay(data.len(), out.len());
    if pairs != t.pairs {
        // The static codes' tables were built for the other choice: with pairs, they could
        // fail where libdeflate's error differs, without the retry below.
        t.static_loaded = false;
        t.pairs = pairs;
    }
    let result = inflate_dispatch(t, data, out);
    if t.pairs && result == Err(DecompressionError::InsufficientSpace) {
        // A pair needs room for both its literals at once, where libdeflate writes the
        // first, then refills, and only then finds no room: on input that has run out, that
        // refill fails first. Decoding again literal by literal gives libdeflate's error.
        t.pairs = false;
        t.static_loaded = false;
        return inflate_dispatch(t, data, out);
    }
    result
}

/// [`inflate_body`], built for the CPU's features.
#[inline(always)]
fn inflate_dispatch(
    t: &mut Tables,
    data: &[u8],
    out: &mut [u8],
) -> Result<Inflated, DecompressionError> {
    #[cfg(target_arch = "x86_64")]
    if !cfg!(feature = "force-portable") && is_x86_feature_detected!("bmi2") {
        // SAFETY: the CPU has BMI2.
        return unsafe { inflate_bmi2(t, data, out) };
    }
    // SAFETY: without BMI2, nothing is required.
    unsafe { inflate_body::<false>(t, data, out) }
}

/// [`inflate_body`] compiled with BMI2, as libdeflate also builds its decompressor:
/// `bzhi` extracts the low bits of the buffer, and `shrx` shifts without tying up `cl`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "bmi2")]
fn inflate_bmi2(
    t: &mut Tables,
    data: &[u8],
    out: &mut [u8],
) -> Result<Inflated, DecompressionError> {
    // SAFETY: this function runs only on a CPU with BMI2.
    unsafe { inflate_body::<true>(t, data, out) }
}

/// The decompressor, inlined into each [`inflate`] variant so it compiles for its target
/// features.
///
/// # Safety
///
/// With `BMI2`, the CPU must have BMI2.
#[inline(always)]
unsafe fn inflate_body<const BMI2: bool>(
    t: &mut Tables,
    data: &[u8],
    out: &mut [u8],
) -> Result<Inflated, DecompressionError> {
    let mut bits = Bits::new(data);
    let mut pos = 0;
    loop {
        let block = start_block(t, &mut bits, out, &mut pos)?;
        if let Some(litlen_mask) = block.huffman
            // SAFETY: as the caller promises.
            && !unsafe { fastloop_for::<BMI2>(t, &mut bits, out, &mut pos, litlen_mask) }?
        {
            generic_loop(t, &mut bits, out, &mut pos, litlen_mask)?;
        }
        if block.is_final {
            break;
        }
    }
    Ok(Inflated {
        consumed: bits.consumed()?,
        written: pos,
    })
}

/// [`fastloop`] for [`inflate_body`]: with BMI2, a function of its own, compiled for it.
/// Inlined, its registers were allocated with the rest of the decoder's, and LLVM kept the
/// table pointer on the stack and reloaded it every symbol, as on Arm: 2% slower on a Xeon
/// Gold 5418Y. Out of line it is within a point of the inlined loop before that refactor.
///
/// # Safety
///
/// With `BMI2`, the CPU must have BMI2.
#[inline(always)]
unsafe fn fastloop_for<const BMI2: bool>(
    t: &Tables,
    bits: &mut Bits<'_>,
    out: &mut [u8],
    pos: &mut usize,
    litlen_mask: u64,
) -> Result<bool, DecompressionError> {
    #[cfg(target_arch = "x86_64")]
    if BMI2 {
        // SAFETY: with `BMI2`, the caller promises the CPU has BMI2.
        return unsafe { fastloop_bmi2(t, bits, out, pos, litlen_mask) };
    }
    fastloop(t, bits, out, pos, litlen_mask)
}

/// [`fastloop`] compiled with BMI2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "bmi2")]
#[inline(never)]
fn fastloop_bmi2(
    t: &Tables,
    bits: &mut Bits<'_>,
    out: &mut [u8],
    pos: &mut usize,
    litlen_mask: u64,
) -> Result<bool, DecompressionError> {
    fastloop(t, bits, out, pos, litlen_mask)
}

/// A block's header, read by [`start_block`].
struct BlockStart {
    is_final: bool,
    /// For a Huffman block, the literal/length table's index mask; `None` for a stored
    /// block, which [`start_block`] has already copied.
    huffman: Option<u64>,
}

/// Reads a block's header: copies a stored block's bytes, or builds a Huffman block's
/// tables.
#[inline(always)]
fn start_block(
    t: &mut Tables,
    bits: &mut Bits<'_>,
    out: &mut [u8],
    pos: &mut usize,
) -> Result<BlockStart, DecompressionError> {
    use DecompressionError::{BadData, InsufficientSpace};

    bits.refill()?;
    let header = bits.buf();
    let is_final = header & 1 != 0;
    match (header >> 1) & 3 {
        0 => {
            // Stored: LEN, NLEN, then LEN bytes, from the next byte boundary.
            bits.consume(3);
            let rest = bits.align_to_byte()?;
            let (len, nlen) = match rest {
                [a, b, c, d, ..] => (u16::from_le_bytes([*a, *b]), u16::from_le_bytes([*c, *d])),
                _ => return Err(BadData),
            };
            if len != !nlen {
                return Err(BadData);
            }
            let len = usize::from(len);
            let dst = out.get_mut(*pos..*pos + len).ok_or(InsufficientSpace)?;
            dst.copy_from_slice(rest.get(4..4 + len).ok_or(BadData)?);
            bits.skip(4 + len);
            *pos += len;
            return Ok(BlockStart {
                is_final,
                huffman: None,
            });
        }
        1 => {
            bits.consume(3);
            if !t.static_loaded {
                t.static_loaded = true;
                t.lens[..144].fill(8);
                t.lens[144..256].fill(9);
                t.lens[256..280].fill(7);
                t.lens[280..288].fill(8);
                t.lens[288..320].fill(5);
                build_tables(t, NUM_LITLEN_SYMS, NUM_OFFSET_SYMS)?;
            }
        }
        2 => {
            let (num_litlen, num_offset) = read_code_lengths(t, bits)?;
            build_tables(t, num_litlen, num_offset)?;
        }
        _ => return Err(BadData),
    }
    Ok(BlockStart {
        is_final,
        huffman: Some(mask(t.litlen_bits)),
    })
}

/// Decodes the rest of a Huffman block with every bounds check, near the ends of the
/// buffers where [`fastloop`] stops.
#[inline(always)]
fn generic_loop(
    t: &Tables,
    bits: &mut Bits<'_>,
    out: &mut [u8],
    pos: &mut usize,
    litlen_mask: u64,
) -> Result<(), DecompressionError> {
    use DecompressionError::{BadData, InsufficientSpace};

    loop {
        // A refill covers a whole symbol: a length (up to 15 + 5 bits) and its offset
        // (up to 15 + 13).
        const _: () = assert!(CONSUMABLE >= 15 + 5 + 15 + 13);
        bits.refill()?;
        let mut entry = t.litlen[(bits.buf() & litlen_mask) as usize];
        let mut saved = bits.buf();
        bits.consume_entry(entry);
        // A pair's second literal can look like the flag, so check it is no literal.
        if entry & (LITERAL | SUBTABLE_POINTER) == SUBTABLE_POINTER {
            entry = t.litlen[subtable_index(entry, bits.buf())];
            saved = bits.buf();
            bits.consume_entry(entry);
        }
        let base = (entry >> 16) as usize;
        if entry & LITERAL != 0 {
            if entry & PAIR != 0 {
                // Where libdeflate would write the first and then run out of room,
                // this fails at once; `inflate` then decodes the stream again
                // without pairs, so that the error is libdeflate's.
                let dst = out.get_mut(*pos..*pos + 2).ok_or(InsufficientSpace)?;
                dst.copy_from_slice(&[base as u8, (entry >> 8) as u8]);
                *pos += 2;
            } else {
                *out.get_mut(*pos).ok_or(InsufficientSpace)? = base as u8;
                *pos += 1;
            }
            continue;
        }
        if entry & END_OF_BLOCK != 0 {
            return Ok(());
        }
        let length = base + extra_bits(entry, saved);
        if length > out.len() - *pos {
            return Err(InsufficientSpace);
        }

        let mut entry = t.offset[(bits.buf() & mask(OFFSET_TABLEBITS)) as usize];
        if entry & EXCEPTIONAL != 0 {
            bits.consume(OFFSET_TABLEBITS);
            entry = t.offset[subtable_index(entry, bits.buf())];
        }
        let offset = (entry >> 16) as usize + extra_bits(entry, bits.buf());
        bits.consume_entry(entry);
        if offset > *pos {
            return Err(BadData);
        }
        copy_match(out, *pos, offset, length);
        *pos += length;
    }
}

/// The most bits an offset takes: a 15-bit codeword and 13 extra bits.
const OFFSET_MAXBITS: u32 = 15 + 13;

/// The most bits an offset takes without a subtable.
const OFFSET_MAXFASTBITS: u32 = OFFSET_TABLEBITS + 13;

/// The bits beyond those counted consumable that the buffer always holds in the fastloop,
/// for looking up the next entry: a fast refill fills all 64 but counts at most 63.
const PRELOAD_SLACK: u32 = 1;

/// The most input one fastloop iteration reads past where it starts: the bytes its refills
/// can count (a full buffer, two literals, a length, and an offset), and the word the last
/// one reads.
const FASTLOOP_MAX_READ: usize =
    (63 + 2 * LITLEN_TABLEBITS as usize + 15 + 5 + 15 + 13).div_ceil(8) + 8;

/// The most output one fastloop iteration writes past where it starts: two literal entries
/// of up to two literals each, the longest match, and the slack its copy may write past it.
const FASTLOOP_MAX_WRITTEN: usize = 4 + 258 + COPY_SLACK;

/// The fast refill and consumption of [`Bits`], in locals the fastloop can keep in
/// registers, with a pointer for the input.
#[derive(Clone, Copy)]
struct FastBits {
    next: *const u8,
    buf: u64,
    left: u32,
}

impl FastBits {
    /// [`Bits::refill`] without its bounds checks.
    ///
    /// # Safety
    ///
    /// The 8 bytes at `next` must be readable.
    #[inline(always)]
    unsafe fn refill(&mut self) {
        // SAFETY: as the caller promises; the read is unaligned.
        let word = u64::from_le(unsafe { self.next.cast::<u64>().read_unaligned() });
        self.buf |= word << (self.left & 63);
        // SAFETY: at most 7 of the 8 readable bytes.
        self.next = unsafe { self.next.add(7 - ((self.left as usize >> 3) & 7)) };
        self.left |= 56;
    }

    #[inline(always)]
    fn consume(&mut self, n: u32) {
        self.buf >>= n & 63;
        self.left = self.left.wrapping_sub(n);
    }

    #[inline(always)]
    fn left(&self) -> u32 {
        u32::from(self.left as u8)
    }
}

/// Decodes literals and matches while one iteration's worst case fits in what is left of
/// both buffers, with no bounds checks on them, as libdeflate's fastloop does; each
/// iteration ([`Fast::step`]) looks up the next entry from bits already in the buffer before
/// its refill and match copy, so they overlap the lookup. Returns whether it reached the end
/// of the block; if not, the generic loop goes on from where it stopped.
///
/// On Arm it is a function of its own: inlined into [`inflate_body`], LLVM kept its table
/// pointer and output bound on the stack and reloaded them every symbol (3-5% slower on an
/// M4 Pro).
#[cfg_attr(target_arch = "aarch64", inline(never))]
#[cfg_attr(not(target_arch = "aarch64"), inline(always))]
fn fastloop(
    t: &Tables,
    bits: &mut Bits<'_>,
    out: &mut [u8],
    pos: &mut usize,
    litlen_mask: u64,
) -> Result<bool, DecompressionError> {
    if !fast_room(bits, out.len(), *pos) {
        return Ok(false);
    }
    // SAFETY: `fast_room` holds, and nothing changes `t` while `f` is used.
    let mut f = unsafe { Fast::new(t, bits, out, *pos, litlen_mask) };
    let mut end_of_block = false;
    while f.room() {
        // SAFETY: `room` holds.
        if let Step::EndOfBlock = unsafe { f.step() }? {
            end_of_block = true;
            break;
        }
    }
    f.save(bits, pos);
    Ok(end_of_block)
}

/// Whether [`fastloop`] can start: one iteration's worst case fits in both buffers.
fn fast_room(bits: &Bits<'_>, out_len: usize, pos: usize) -> bool {
    let (data, in_pos, _, _) = bits.fast_state();
    data.len() - in_pos > FASTLOOP_MAX_READ && pos < out_len.saturating_sub(FASTLOOP_MAX_WRITTEN)
}

/// One stream's fastloop state, in locals that stay in registers: [`fastloop`]'s, and each
/// of the batch decoder's two lanes'.
struct Fast {
    /// The stream's tables, which nothing changes while the fastloop runs.
    t: *const Tables,
    litlen_mask: u64,
    b: FastBits,
    /// The next literal/length entry, looked up but not consumed.
    entry: u32,
    out_next: *mut u8,
    out_base: *mut u8,
    in_fast_end: *const u8,
    out_fast_end: *mut u8,
}

/// What a [`Fast::step`] reached.
enum Step {
    Go,
    EndOfBlock,
}

impl Fast {
    /// The state at `bits` and `pos` in `out`, with the first entry looked up.
    ///
    /// # Safety
    ///
    /// [`fast_room`] must hold for `bits`, `out` and `pos`, and nothing may change `t` while
    /// the result is used.
    #[inline(always)]
    unsafe fn new(
        t: &Tables,
        bits: &Bits<'_>,
        out: &mut [u8],
        pos: usize,
        litlen_mask: u64,
    ) -> Self {
        let (data, in_pos, buf, left) = bits.fast_state();
        debug_assert!(litlen_mask < 1 << LITLEN_TABLEBITS, "mask {litlen_mask:#x}");
        // SAFETY: `fast_room` puts more than `FASTLOOP_MAX_READ` bytes of `data` after
        // `in_pos` and more than `FASTLOOP_MAX_WRITTEN` of `out` after `pos`.
        unsafe {
            let out_base = out.as_mut_ptr();
            let mut f = Self {
                t,
                litlen_mask,
                b: FastBits {
                    next: data.as_ptr().add(in_pos),
                    buf,
                    left,
                },
                entry: 0,
                out_next: out_base.add(pos),
                out_base,
                in_fast_end: data.as_ptr().add(data.len() - FASTLOOP_MAX_READ),
                out_fast_end: out_base.add(out.len() - FASTLOOP_MAX_WRITTEN),
            };
            f.b.refill();
            f.entry = f.lookup(f.b.buf);
            f
        }
    }

    /// Whether another [`Self::step`] can run.
    #[inline(always)]
    fn room(&self) -> bool {
        self.b.next < self.in_fast_end && self.out_next < self.out_fast_end
    }

    #[inline(always)]
    fn lookup(&self, buf: u64) -> u32 {
        const _: () = assert!(1 << LITLEN_TABLEBITS <= LITLEN_ENOUGH);
        // SAFETY: `t` is live and unchanged while the fastloop runs; the main table has
        // `1 << t.litlen_bits` entries, `litlen_mask` selects that many bits, and
        // `litlen_bits <= LITLEN_TABLEBITS`, whose main table fits in `LITLEN_ENOUGH`.
        unsafe {
            *(*self.t)
                .litlen
                .get_unchecked((buf & self.litlen_mask) as usize)
        }
    }

    /// One iteration of the fastloop: a literal entry or up to three, or a match.
    ///
    /// # Safety
    ///
    /// [`Self::room`] must hold.
    #[inline(always)]
    unsafe fn step(&mut self) -> Result<Step, DecompressionError> {
        // SAFETY: `t` is live and unchanged while the fastloop runs.
        let t = unsafe { &*self.t };
        let b = &mut self.b;
        let mut entry = self.entry;
        // SAFETY: `room` holds: `next` is more than `FASTLOOP_MAX_READ` bytes from the end
        // of the input and `out_next` more than `FASTLOOP_MAX_WRITTEN` from the end of the
        // output, and those are what one iteration can read and write at most: its refills
        // advance `next` by the bits it consumes and read a word past it, and it writes at
        // most two literals and a match with its copy's slack. Matches start at most
        // `out_next - out_base` back, as checked.
        unsafe {
            let mut saved = b.buf;
            b.consume(entry);
            if entry & LITERAL != 0 {
                // Up to two more literals before the next refill, while there is room in
                // the bits for them and whatever follows. Each may be a pair.
                let lit = entry;
                entry = *t.litlen.get_unchecked((b.buf & self.litlen_mask) as usize);
                saved = b.buf;
                b.consume(entry);
                self.out_next = push_literal(self.out_next, lit);
                if entry & LITERAL != 0 {
                    let lit = entry;
                    entry = *t.litlen.get_unchecked((b.buf & self.litlen_mask) as usize);
                    saved = b.buf;
                    b.consume(entry);
                    self.out_next = push_literal(self.out_next, lit);
                    if entry & LITERAL != 0 {
                        // A third, in place of this iteration's length.
                        let lit = entry;
                        self.entry = *t.litlen.get_unchecked((b.buf & self.litlen_mask) as usize);
                        b.refill();
                        self.out_next = push_literal(self.out_next, lit);
                        return Ok(Step::Go);
                    }
                }
            }
            if entry & EXCEPTIONAL != 0 {
                if entry & END_OF_BLOCK != 0 {
                    return Ok(Step::EndOfBlock);
                }
                entry = t.litlen[subtable_index(entry, b.buf)];
                saved = b.buf;
                b.consume(entry);
                if entry & LITERAL != 0 {
                    let lit = (entry >> 16) as u8;
                    self.entry = *t.litlen.get_unchecked((b.buf & self.litlen_mask) as usize);
                    b.refill();
                    *self.out_next = lit;
                    self.out_next = self.out_next.add(1);
                    return Ok(Step::Go);
                }
                if entry & END_OF_BLOCK != 0 {
                    return Ok(Step::EndOfBlock);
                }
            }
            let length = (entry >> 16) as usize + extra_bits(entry, saved);

            // The offset, refilling first if it and the next lookup might not fit.
            entry = t.offset[(b.buf & mask(OFFSET_TABLEBITS)) as usize];
            if entry & EXCEPTIONAL != 0 {
                if b.left() < OFFSET_MAXBITS + LITLEN_TABLEBITS - PRELOAD_SLACK {
                    b.refill();
                }
                b.consume(OFFSET_TABLEBITS);
                entry = t.offset[subtable_index(entry, b.buf)];
            } else if b.left() < OFFSET_MAXFASTBITS + LITLEN_TABLEBITS - PRELOAD_SLACK {
                b.refill();
            }
            let saved = b.buf;
            b.consume(entry);
            let offset = (entry >> 16) as usize + extra_bits(entry, saved);
            if offset > self.out_next.offset_from_unsigned(self.out_base) {
                return Err(DecompressionError::BadData);
            }

            self.entry = *t.litlen.get_unchecked((b.buf & self.litlen_mask) as usize);
            b.refill();
            copy_words(self.out_next, offset, length);
            self.out_next = self.out_next.add(length);
            Ok(Step::Go)
        }
    }

    /// Writes the place reached back to `bits` and `pos`.
    #[inline(always)]
    fn save(self, bits: &mut Bits<'_>, pos: &mut usize) {
        let (data, ..) = bits.fast_state();
        // SAFETY: both pointers stayed within their buffers, which `room` keeps them in.
        unsafe {
            bits.set_fast_state(
                self.b.next.offset_from_unsigned(data.as_ptr()),
                self.b.buf,
                self.b.left,
            );
            *pos = self.out_next.offset_from_unsigned(self.out_base);
        }
    }
}

/// The bytes a match copy works in.
const WORD: usize = 8;

/// The bytes a match copy may write past the match's end: it copies five words before
/// checking whether it is done, as libdeflate does, since most matches are that short and
/// a loop exit that depends on the length predicts badly.
const COPY_SLACK: usize = 5 * WORD - 1;

/// Copies the `length` bytes `offset` back from `pos` in `out` to `pos`. The copy may
/// overlap itself: then it repeats the last `offset` bytes.
///
/// With room for [`COPY_SLACK`] bytes past the end, it copies a word at a time, as
/// libdeflate's fastloop does, writing past the match bytes that later output overwrites.
///
/// # Panics
///
/// If `offset` is 0 or more than `pos`, or the match does not fit in `out`.
#[inline(always)]
fn copy_match(out: &mut [u8], pos: usize, offset: usize, length: usize) {
    assert!(0 < offset && offset <= pos, "offset {offset} at {pos}");
    let end = pos + length;
    if end + COPY_SLACK > out.len() {
        for i in pos..end {
            out[i] = out[i - offset];
        }
        return;
    }
    // SAFETY: checked above: `offset <= pos`, and `COPY_SLACK` bytes of room past `end`.
    unsafe { copy_words(out.as_mut_ptr().add(pos), offset, length) };
}

/// The word-at-a-time part of [`copy_match`]: copies the `length` bytes `offset` back from
/// `dst` to `dst`, writing up to [`COPY_SLACK`] bytes past them.
///
/// # Safety
///
/// `dst - offset..dst + length + COPY_SLACK` must be valid for reads and writes, with
/// `offset > 0` and `length >= 3`.
#[inline(always)]
unsafe fn copy_words(mut dst: *mut u8, offset: usize, length: usize) {
    // SAFETY: every access is in `dst - offset..dst + length + COPY_SLACK`, as the caller
    // promises. The stores before the first check of `dst < end` reach at most
    // `dst + 5 * WORD <= end + COPY_SLACK` (a match is at least 3 bytes), and each later one
    // starts below `end`, so ends by `end + WORD - 1`. Each load is `offset` or more before
    // its store, so starts at `dst - offset` or later.
    unsafe {
        let end = dst.add(length);
        let mut src = dst.sub(offset);
        let copy_word = |src: *const u8, dst: *mut u8| {
            dst.cast::<u64>()
                .write_unaligned(src.cast::<u64>().read_unaligned());
        };
        if offset >= WORD {
            // Each word's source is already written.
            loop {
                for i in 0..5 {
                    copy_word(src.add(i * WORD), dst.add(i * WORD));
                }
                src = src.add(5 * WORD);
                dst = dst.add(5 * WORD);
                if dst >= end {
                    break;
                }
            }
        } else if offset == 1 {
            // A run of one byte, four words at a time (which compiles to two 16-byte stores).
            let word = u64::from(*src) * 0x0101_0101_0101_0101;
            loop {
                for i in 0..4 {
                    dst.add(i * WORD).cast::<u64>().write_unaligned(word);
                }
                dst = dst.add(4 * WORD);
                if dst >= end {
                    break;
                }
            }
        } else {
            // Each word's first `offset` bytes are already written, which is as far as it
            // advances; the next word rewrites the rest. A match is at least 3 bytes, and
            // `offset >= 2`, so two words always help.
            copy_word(src, dst);
            copy_word(src.add(offset), dst.add(offset));
            src = src.add(2 * offset);
            dst = dst.add(2 * offset);
            while dst < end {
                copy_word(src, dst);
                copy_word(src.add(offset), dst.add(offset));
                src = src.add(2 * offset);
                dst = dst.add(2 * offset);
            }
        }
    }
}

/// Writes a literal entry's literal, or both of a [`PAIR`]'s, at `dst`, and returns where
/// the output goes on.
///
/// # Safety
///
/// Two bytes at `dst` must be writable: it writes two either way, without a branch.
#[inline(always)]
unsafe fn push_literal(dst: *mut u8, entry: u32) -> *mut u8 {
    let two = ((entry >> 16) & 0xff) | (entry & 0xff00);
    // SAFETY: as the caller promises.
    unsafe {
        dst.cast::<u16>().write_unaligned((two as u16).to_le());
        dst.add(1 + ((entry >> 30) & 1) as usize)
    }
}

/// Whether literal pairs are worth building for a stream of `compressed` bytes that
/// inflates to about `inflated`.
///
/// Building pairs costs a few thousand cycles per table and saves a lookup each time one
/// is decoded, so it pays where each table decodes many symbols: data that compresses
/// poorly, like long reads' qualities and tags, rather than highly compressible data like
/// BCF and VCF, whose tables each cover few symbols. Compressed bits per output byte is
/// known before decoding; the cut-off sits between the two kinds (about 2-5:1 against
/// 12:1 and more). `inflated` is the output buffer's size, exact for BGZF; a larger buffer
/// only makes pairs less likely.
fn pairs_pay(compressed: usize, inflated: usize) -> bool {
    compressed.saturating_mul(PAIRS_MAX_RATIO) >= inflated
}

/// The compression ratio above which literal pairs are not built. Arm's cores build tables
/// faster relative to decoding, so pairs pay to a higher ratio there: 16 was 3-4% faster
/// than 8 on Illumina BAM on an M4 Pro, and no slower on long reads or VCF.
const PAIRS_MAX_RATIO: usize = if cfg!(target_arch = "aarch64") { 16 } else { 8 };

/// The index of the subtable entry that the subtable pointer `entry` leads to, for the bits
/// in `buf` after the main table's.
fn subtable_index(entry: u32, buf: u64) -> usize {
    (entry >> 16) as usize + (buf & mask((entry >> 8) & 0x3f)) as usize
}

/// The value of a length's or offset's extra bits, which follow its codeword in `buf`.
fn extra_bits(entry: u32, buf: u64) -> usize {
    ((buf & mask(entry & 0xff)) >> ((entry >> 8) & 0xff)) as usize
}

/// Reads a dynamic block's header: the precode, then the literal/length and offset code
/// lengths into `t.lens`. Returns how many of each there are.
fn read_code_lengths(
    t: &mut Tables,
    bits: &mut Bits<'_>,
) -> Result<(usize, usize), DecompressionError> {
    #[cfg(target_arch = "x86_64")]
    if !cfg!(feature = "force-portable") && is_x86_feature_detected!("bmi2") {
        // SAFETY: the CPU has BMI2.
        return unsafe { read_code_lengths_bmi2(t, bits) };
    }
    read_code_lengths_body(t, bits)
}

/// [`read_code_lengths_body`] with BMI2: its loop is a chain of lookups and shifts, and
/// `shrx` shortens the shifts. Chosen here rather than by the caller, so that the
/// decoder's own build is unchanged.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "bmi2")]
fn read_code_lengths_bmi2(
    t: &mut Tables,
    bits: &mut Bits<'_>,
) -> Result<(usize, usize), DecompressionError> {
    read_code_lengths_body(t, bits)
}

#[inline(always)]
fn read_code_lengths_body(
    t: &mut Tables,
    bits: &mut Bits<'_>,
) -> Result<(usize, usize), DecompressionError> {
    use DecompressionError::BadData;

    // A local copy stays in registers; through the reference, a possible panic makes the
    // compiler store it back on every symbol.
    let mut b = bits.clone();

    // BFINAL and BTYPE, HLIT, HDIST, HCLEN, and the first precode length.
    const _: () = assert!(CONSUMABLE >= 3 + 5 + 5 + 4 + 3);
    let header = b.buf();
    let num_litlen = 257 + ((header >> 3) & 0x1f) as usize;
    let num_offset = 1 + ((header >> 8) & 0x1f) as usize;
    let num_precode = 4 + ((header >> 13) & 0xf) as usize;
    t.static_loaded = false;

    // A refill leaves room for the other 18 precode lengths.
    const _: () = assert!(CONSUMABLE >= 3 * (NUM_PRECODE_SYMS as u32 - 1));
    t.precode_lens[PRECODE_ORDER[0]] = ((header >> 17) & 7) as u8;
    b.consume(20);
    b.refill()?;
    for &sym in &PRECODE_ORDER[1..num_precode] {
        t.precode_lens[sym] = b.peek(3) as u8;
        b.consume(3);
    }
    for &sym in &PRECODE_ORDER[num_precode..] {
        t.precode_lens[sym] = 0;
    }
    table::build(
        &mut t.precode,
        &t.precode_lens,
        &PRECODE_RESULTS,
        TableBits::Fixed(PRECODE_TABLEBITS),
        MAX_PRE_CODEWORD_LEN,
        &mut t.sorted,
        Pairs::None,
    )
    .ok_or(BadData)?;

    // The code lengths, three to a refill: a refill leaves room for three symbols with
    // their repeat counts, so no symbol needs to test whether to refill first, a test that
    // comes out as good as at random.
    const _: () = assert!(CONSUMABLE >= 3 * (MAX_PRE_CODEWORD_LEN + 7));
    let total = num_litlen + num_offset;
    let mut i = 0;
    'lengths: loop {
        b.refill()?;
        for _ in 0..3 {
            if i >= total {
                break 'lengths;
            }
            let entry = t.precode[b.peek(MAX_PRE_CODEWORD_LEN) as usize];
            b.consume_entry(entry);
            let presym = (entry >> 16) as u8;
            let (value, count) = match presym {
                0..16 => {
                    t.lens[i] = presym;
                    i += 1;
                    continue;
                }
                16 => {
                    // Repeat the previous length 3-6 times.
                    let prev = *t.lens[..i].last().ok_or(BadData)?;
                    let count = 3 + b.peek(2) as usize;
                    b.consume(2);
                    (prev, count)
                }
                17 => {
                    let count = 3 + b.peek(3) as usize;
                    b.consume(3);
                    (0, count)
                }
                _ => {
                    let count = 11 + b.peek(7) as usize;
                    b.consume(7);
                    (0, count)
                }
            };
            t.lens[i..i + count].fill(value);
            i += count;
        }
    }
    // A run must not go past the last length (as libdeflate and zlib check).
    if i != total {
        return Err(BadData);
    }
    *bits = b;
    Ok((num_litlen, num_offset))
}

/// Builds the offset and literal/length tables from `t.lens`.
///
/// Never inlined: inlined into a decoder at either call site, it slows the hot loop by ~5%
/// (code layout and register allocation), more than the call costs.
#[inline(never)]
fn build_tables(
    t: &mut Tables,
    num_litlen: usize,
    num_offset: usize,
) -> Result<(), DecompressionError> {
    // A wider table where it holds pairs, so that more of them fit.
    let bits = if t.pairs {
        LITLEN_TABLEBITS
    } else {
        LITLEN_TABLEBITS_PLAIN
    };
    let (litlen_lens, rest) = t.lens.split_at(num_litlen);
    table::build(
        &mut t.offset,
        &rest[..num_offset],
        &OFFSET_RESULTS,
        TableBits::Fixed(OFFSET_TABLEBITS),
        MAX_CODEWORD_LEN,
        &mut t.sorted,
        Pairs::None,
    )
    .ok_or(DecompressionError::BadData)?;
    t.litlen_bits = table::build(
        &mut t.litlen,
        litlen_lens,
        &LITLEN_RESULTS,
        TableBits::AtMost(bits),
        MAX_CODEWORD_LEN,
        &mut t.sorted,
        if t.pairs { Pairs::Add } else { Pairs::None },
    )
    .ok_or(DecompressionError::BadData)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The plain build must agree with the one [`inflate`] picks, which the differential
    /// tests check against libdeflate; on a CPU with BMI2 they are different code.
    #[test]
    fn plain_build_matches_dispatched() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let (mut a, mut b) = (Tables::try_new().unwrap(), Tables::try_new().unwrap());
        for level in [0, 1, 6, 12] {
            let mut c =
                libdeflater::Compressor::new(libdeflater::CompressionLvl::new(level).unwrap());
            for len in [0, 1, 100, 4096, 65_536] {
                let data: Vec<u8> = (0..len)
                    .map(|i| b"ACGT\tchr1"[(next() as usize + i) % 9])
                    .collect();
                let mut z = vec![0; c.deflate_compress_bound(len)];
                let n = c.deflate_compress(&data, &mut z).unwrap();
                z.truncate(n);
                for trial in 0..300 {
                    let mut input = z.clone();
                    if trial > 0 && !input.is_empty() {
                        let i = next() as usize % input.len();
                        input[i] ^= 1 << (next() % 8);
                    }
                    let (mut x, mut y) = (vec![0; len + 64], vec![0; len + 64]);
                    // SAFETY: without BMI2, nothing is required.
                    let rx = unsafe { inflate_body::<false>(&mut a, &input, &mut x) };
                    let ry = inflate(&mut b, &input, &mut y);
                    assert_eq!(rx, ry, "level {level} len {len} trial {trial}");
                    if let Ok(r) = rx {
                        assert!(
                            x[..r.written] == y[..r.written],
                            "level {level} len {len} trial {trial}"
                        );
                    }
                }
            }
        }
    }
}
