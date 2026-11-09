//! Adler-32 without SIMD: libdeflate's `adler32_generic` and `ADLER32_CHUNK` (MIT), which add
//! four bytes a step into separate sums, for more instruction-level parallelism than adding one
//! byte at a time.

use super::{DIVISOR, MAX_CHUNK_LEN};

/// Adler-32 four bytes at a time.
pub(super) fn generic(adler: u32, data: &[u8]) -> u32 {
    let (mut s1, mut s2) = (adler & 0xffff, adler >> 16);
    for c in data.chunks(MAX_CHUNK_LEN & !3) {
        chunk(&mut s1, &mut s2, c);
    }
    (s2 << 16) | s1
}

/// Adds `data` to `s1` and `s2`, then reduces them modulo [`DIVISOR`].
///
/// The sums must not overflow: `data`, with what was added since they were last reduced, must
/// be at most [`MAX_CHUNK_LEN`] bytes.
#[inline(always)]
pub(super) fn chunk(s1: &mut u32, s2: &mut u32, data: &[u8]) {
    let (quads, rest) = data.as_chunks::<4>();
    if !quads.is_empty() {
        // Each step adds the `s1` it starts with four times to `s2`, and its bytes 4, 3, 2
        // and 1 times: sum those apart and weight them at the end.
        let (mut s1_sum, mut b0, mut b1, mut b2, mut b3) = (0, 0, 0, 0, 0);
        for q in quads {
            let [x0, x1, x2, x3] = q.map(u32::from);
            s1_sum += *s1;
            *s1 += x0 + x1 + x2 + x3;
            b0 += x0;
            b1 += x1;
            b2 += x2;
            b3 += x3;
        }
        *s2 += 4 * (s1_sum + b0) + 3 * b1 + 2 * b2 + b3;
    }
    for &b in rest {
        *s1 += u32::from(b);
        *s2 += *s1;
    }
    *s1 %= DIVISOR;
    *s2 %= DIVISOR;
}
