//! Writing DEFLATE's bitstream: bits go into a 64-bit buffer, low bits first, and whole bytes
//! go out eight at a time.

/// The output and the bits not yet written to it. A write that does not fit sets `overflow`
/// and writes nothing more; the compression then fails with no partial result.
pub(super) struct BitWriter<'a> {
    out: &'a mut [u8],
    pos: usize,
    /// `out.len() - 7`, or 0: below it, 8 bytes fit from `pos`. One compare against it is
    /// the whole check for the fast flush.
    fast_end: usize,
    bitbuf: u64,
    /// Bits in `bitbuf`: at most 7 after a flush, and at most 63 ever.
    bitcount: u32,
    overflow: bool,
}

impl<'a> BitWriter<'a> {
    pub(super) fn new(out: &'a mut [u8]) -> Self {
        Self {
            fast_end: out.len().saturating_sub(7),
            out,
            pos: 0,
            bitbuf: 0,
            bitcount: 0,
            overflow: false,
        }
    }

    /// Adds the low `n` bits of `bits`, whose other bits are 0. The caller flushes often
    /// enough that the buffer never holds more than 63 bits.
    #[inline(always)]
    pub(super) fn add(&mut self, bits: u32, n: u32) {
        debug_assert!(n <= 32 && (n == 32 || bits >> n == 0) && self.bitcount + n <= 63);
        self.bitbuf |= u64::from(bits) << self.bitcount;
        self.bitcount += n;
    }

    /// Writes the buffer's whole bytes.
    #[inline(always)]
    pub(super) fn flush(&mut self) {
        if self.pos < self.fast_end {
            // SAFETY: `pos < out.len() - 7`, so 8 bytes from `pos` are in `out`.
            unsafe {
                self.out
                    .as_mut_ptr()
                    .add(self.pos)
                    .cast::<u64>()
                    .write_unaligned(self.bitbuf.to_le());
            }
            self.pos += (self.bitcount >> 3) as usize;
            self.bitbuf >>= self.bitcount & !7;
            self.bitcount &= 7;
        } else {
            crate::cold_path();
            // The slow path gets the state by value: were it handed `self`, the compiler
            // would keep the buffer in memory, not registers, in every loop that flushes.
            let (pos, bitbuf, bitcount, fits) =
                flush_slow(self.out, self.pos, self.bitbuf, self.bitcount);
            self.pos = pos;
            self.bitbuf = bitbuf;
            self.bitcount = bitcount;
            self.overflow |= !fits;
        }
    }

    /// Bits buffered, at most 7 between blocks.
    pub(super) fn bitcount(&self) -> u32 {
        self.bitcount
    }

    /// Bytes of output left.
    pub(super) fn room(&self) -> usize {
        self.out.len() - self.pos
    }

    /// Bits written so far, buffered ones included.
    pub(super) fn bits_written(&self) -> usize {
        8 * self.pos + self.bitcount as usize
    }

    pub(super) fn set_overflow(&mut self) {
        self.overflow = true;
    }

    pub(super) fn overflowed(&self) -> bool {
        self.overflow
    }

    /// Writes a stored block's header (BFINAL, then BTYPE 0), pads to a byte, and writes
    /// LEN, NLEN, and `data`, at most 65,535 bytes.
    pub(super) fn write_stored(&mut self, is_final: bool, data: &[u8]) {
        debug_assert!(data.len() <= usize::from(u16::MAX) && self.bitcount <= 7);
        // After an overflow, nothing more is written, and the buffered bits stay as they
        // were: a header added for every block that did not fit would overrun them.
        if self.overflow {
            return;
        }
        let bits = self.bitbuf | u64::from(is_final) << self.bitcount;
        let header_bytes = (self.bitcount + 3).div_ceil(8) as usize;
        let len = data.len() as u16;
        let Some(dst) = self
            .out
            .get_mut(self.pos..self.pos + header_bytes + 4 + data.len())
        else {
            self.overflow = true;
            return;
        };
        let (header, rest) = dst.split_at_mut(header_bytes);
        header.copy_from_slice(&bits.to_le_bytes()[..header_bytes]);
        rest[..2].copy_from_slice(&len.to_le_bytes());
        rest[2..4].copy_from_slice(&(!len).to_le_bytes());
        rest[4..].copy_from_slice(data);
        self.pos += header_bytes + 4 + data.len();
        self.bitbuf = 0;
        self.bitcount = 0;
    }

    /// Writes the last partial byte, and returns the length written, or `None` if the output
    /// was too small.
    pub(super) fn finish(self) -> Option<usize> {
        debug_assert!(self.bitcount <= 7);
        let (pos, _, _, fits) = flush_slow(
            self.out,
            self.pos,
            self.bitbuf,
            self.bitcount.next_multiple_of(8),
        );
        (fits && !self.overflow).then_some(pos)
    }
}

/// [`BitWriter::flush`] within 8 bytes of the end, a byte at a time: returns the new
/// position and buffer, and whether everything fit.
#[cold]
#[inline(never)]
fn flush_slow(
    out: &mut [u8],
    mut pos: usize,
    mut bitbuf: u64,
    mut bitcount: u32,
) -> (usize, u64, u32, bool) {
    while bitcount >= 8 {
        let Some(byte) = out.get_mut(pos) else {
            // Nothing more will be kept; drop the bits so the buffer cannot overfill.
            return (pos, 0, 0, false);
        };
        *byte = bitbuf as u8;
        pos += 1;
        bitbuf >>= 8;
        bitcount -= 8;
    }
    (pos, bitbuf, bitcount, true)
}
