//! The input bitstream, low bit first, as libdeflate's `deflate_decompress.c` (MIT) reads
//! it: a 64-bit buffer refilled a word at a time while the input lasts.

use crate::DecompressionError;

/// The most bits `left` counts. One short of 64 lets a refill shift by `left` and mark the
/// whole bytes it added with an OR.
const MAX_LEFT: u32 = 63;

/// The fewest bits a refill leaves to consume: only whole bytes are added.
pub(super) const CONSUMABLE: u32 = MAX_LEFT - 7;

/// The bitstream. Past the end of the input it reads zero bytes, so decoding never has to
/// tell its own lookahead from a real overread; it counts them, and rejects the input if it
/// would need more than a refill's worth, or once a block ends, if it consumed any.
#[derive(Clone, Debug)]
pub(super) struct Bits<'d> {
    data: &'d [u8],
    /// The next byte not yet in `buf`.
    pos: usize,
    /// Bits read but not consumed, the next in bit 0. Bits above `left` may hold input read
    /// ahead, which a later refill ORs in again at the same place.
    buf: u64,
    /// The bits in `buf` that can be consumed, in the low byte. The bits above it may be
    /// garbage, as in libdeflate: consuming a table entry subtracts the whole entry, whose
    /// low byte is the bits it takes, which saves extracting that byte on every symbol.
    left: u32,
    /// The implicit zero bytes read past the end.
    overread: usize,
}

impl<'d> Bits<'d> {
    pub(super) fn new(data: &'d [u8]) -> Self {
        Self {
            data,
            pos: 0,
            buf: 0,
            left: 0,
            overread: 0,
        }
    }

    /// Fills the buffer to at least [`CONSUMABLE`] bits.
    pub(super) fn refill(&mut self) -> Result<(), DecompressionError> {
        if let Some(&word) = self.data[self.pos..].first_chunk::<8>() {
            self.add_word(word);
            return Ok(());
        }
        while self.left() < CONSUMABLE {
            if let Some(&b) = self.data.get(self.pos) {
                self.buf |= u64::from(b) << self.left();
                self.pos += 1;
            } else {
                self.overread += 1;
                if self.overread > 8 {
                    return Err(DecompressionError::BadData);
                }
            }
            self.left = self.left.wrapping_add(8);
        }
        Ok(())
    }

    /// Branchless: takes the word, and counts the whole bytes that fit.
    fn add_word(&mut self, word: [u8; 8]) {
        self.buf |= u64::from_le_bytes(word) << (self.left & 63);
        self.pos += 7 - ((self.left as usize >> 3) & 7);
        self.left |= MAX_LEFT & !7;
    }

    /// The input, the next byte not yet in the buffer, the buffer, and the count of bits in
    /// it, for the fastloop to work on in its own registers. Implicit zero bytes are read
    /// only once the input is used up, so with any input left there are none.
    pub(super) fn fast_state(&self) -> (&'d [u8], usize, u64, u32) {
        (self.data, self.pos, self.buf, self.left)
    }

    /// Takes back the state the fastloop worked on: the next input byte not yet in the
    /// buffer, at most the input's length, the buffer, and the count of bits in it.
    pub(super) fn set_fast_state(&mut self, pos: usize, buf: u64, left: u32) {
        assert!(
            pos <= self.data.len(),
            "input position {pos} of {}",
            self.data.len()
        );
        self.pos = pos;
        self.buf = buf;
        self.left = left;
    }

    /// The next `n` bits, `n < 64`, without consuming them.
    pub(super) fn peek(&self, n: u32) -> u64 {
        self.buf & mask(n)
    }

    /// The whole buffer, for decoding from the bits of an entry as libdeflate does.
    pub(super) fn buf(&self) -> u64 {
        self.buf
    }

    pub(super) fn left(&self) -> u32 {
        u32::from(self.left as u8)
    }

    /// Consumes `n` bits; at most [`left`](Self::left).
    pub(super) fn consume(&mut self, n: u32) {
        debug_assert!(n <= self.left(), "consuming {n} of {} bits", self.left());
        self.buf >>= n;
        self.left = self.left.wrapping_sub(n);
    }

    /// Consumes the bits a table entry takes: its low byte, at most [`left`](Self::left).
    pub(super) fn consume_entry(&mut self, entry: u32) {
        debug_assert!(
            entry & 0xff <= self.left(),
            "consuming {} of {} bits",
            entry & 0xff,
            self.left()
        );
        // Shifts take their count mod 64, and the entry's low byte is under 64; subtracting
        // the whole entry leaves the low byte of `left` exact, with no borrow out of it.
        self.buf >>= entry & 63;
        self.left = self.left.wrapping_sub(entry);
    }

    /// Drops the bits up to the next byte boundary and returns the input from there, for a
    /// stored block, which is byte-aligned. Fails if that would go back past real input.
    pub(super) fn align_to_byte(&mut self) -> Result<&'d [u8], DecompressionError> {
        let unread = (self.left() >> 3) as usize;
        if self.overread > unread {
            return Err(DecompressionError::BadData);
        }
        self.pos -= unread - self.overread;
        self.overread = 0;
        self.buf = 0;
        self.left = 0;
        Ok(&self.data[self.pos..])
    }

    /// Moves past `n` bytes of the input returned by
    /// [`align_to_byte`](Self::align_to_byte).
    pub(super) fn skip(&mut self, n: usize) {
        self.pos += n;
    }

    /// The input bytes consumed, at the end of the stream; fails if that includes implicit
    /// zero bytes.
    pub(super) fn consumed(&self) -> Result<usize, DecompressionError> {
        let unread = (self.left() >> 3) as usize;
        if self.overread > unread {
            return Err(DecompressionError::BadData);
        }
        Ok(self.pos - (unread - self.overread))
    }
}

/// The low `n` bits set, `n < 64`.
pub(super) fn mask(n: u32) -> u64 {
    (1 << n) - 1
}
