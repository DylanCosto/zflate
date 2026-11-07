//! CRC32 with the Armv8 CRC instructions, which 64-bit Arm CPUs that run desktop and server
//! operating systems have (Apple's, AWS Graviton, Ampere): eight bytes an instruction.
//!
//! One instruction waits for the last (about 3 cycles), so the fast version runs three
//! streams at once, over three thirds of each chunk, and joins them with carry-less
//! multiplies (PMULL), as libdeflate's `crc32_arm_crc_pmullcombine` does.
//!
//! Apple's cores have one CRC unit but several PMULL pipes, so there folding the data with
//! PMULL, as the x86 code does, is faster still (about 68 GB/s against 25 on an M1, by
//! libdeflate's measurements): twelve 128-bit accumulators, as libdeflate's
//! `arm/crc32_pmull_wide.h` (MIT) does, with EOR3 from the SHA3 extension. Faster again is
//! to keep the CRC unit busy too, on its own share of the data, as [`crc32_arm_hybrid`]
//! does: about 89 GB/s on an M4 Pro at 64 KiB, where libdeflate's PMULL code makes 76.

use std::arch::aarch64::*;

use super::{Impl, x_pow_mod_g};

/// The implementations this CPU supports, fastest first.
pub(super) fn available() -> Vec<(&'static str, Impl)> {
    let mut impls = Vec::new();
    let mut pmull = None;
    if std::arch::is_aarch64_feature_detected!("crc")
        && std::arch::is_aarch64_feature_detected!("aes")
        && std::arch::is_aarch64_feature_detected!("sha3")
    {
        // SAFETY: returned only once the CRC, PMULL and EOR3 instructions are detected.
        let f: Impl = |crc, data| unsafe { crc32_arm_pmull(crc, data) };
        pmull = Some(("arm-pmull-12way", f));
    }
    // Only Apple's cores are known to fold faster than the CRC unit runs; libdeflate also
    // prefers PMULL on Neoverse V cores, which it tells apart by MIDR_EL1 on Linux.
    if cfg!(target_vendor = "apple") {
        if pmull.is_some() {
            // SAFETY: as for `pmull`: the hybrid needs the same instructions. Short inputs
            // go to the PMULL code from here, a call sooner than through the hybrid.
            let f: Impl = |crc, data| unsafe {
                if data.len() < HYBRID_CHUNK {
                    crc32_arm_pmull(crc, data)
                } else {
                    crc32_arm_hybrid(crc, data)
                }
            };
            impls.push(("arm-hybrid", f));
        }
        impls.extend(pmull.take());
    }
    if std::arch::is_aarch64_feature_detected!("crc")
        && std::arch::is_aarch64_feature_detected!("aes")
    {
        // SAFETY: returned only once the CRC and PMULL instructions are detected.
        let f: Impl = |crc, data| unsafe { crc32_arm_3way(crc, data) };
        impls.push(("arm-crc-3way", f));
    }
    impls.extend(pmull);
    if std::arch::is_aarch64_feature_detected!("crc") {
        // SAFETY: returned only once the CRC instructions are detected.
        let f: Impl = |crc, data| unsafe { crc32_arm(crc, data) };
        impls.push(("arm-crc", f));
    }
    impls
}

/// The CRC instructions work on the register (the inverted CRC), as [`Impl`] does: eight
/// bytes at a time, then the last few one at a time.
#[target_feature(enable = "crc")]
fn crc32_arm(mut crc: u32, data: &[u8]) -> u32 {
    let (words, tail) = data.as_chunks::<8>();
    for word in words {
        crc = __crc32d(crc, u64::from_le_bytes(*word));
    }
    for &byte in tail {
        crc = __crc32b(crc, byte);
    }
    crc
}

/// Each of the three streams' share of a chunk, in bytes.
const STRIDE: usize = 512;

/// The multipliers that move a register `2 * STRIDE` and `STRIDE` bytes on: `x^(8n - 33)
/// mod G` for `n` bytes, the 33 undone by the multiply's one-bit shift and the CRC
/// instruction's 32.
const SHIFT_2: u32 = x_pow_mod_g(8 * 2 * STRIDE as u32 - 33);
const SHIFT_1: u32 = x_pow_mod_g(8 * STRIDE as u32 - 33);

/// The register `crc` moved on by the bytes `multiplier` stands for: the carry-less product,
/// reduced by a CRC instruction.
#[target_feature(enable = "crc,aes")]
fn shift(crc: u32, multiplier: u32) -> u32 {
    let product = vmull_p64(u64::from(crc), u64::from(multiplier));
    __crc32d(0, product as u64)
}

/// Chunks of `3 * STRIDE` bytes as three streams, the first carrying the register in and
/// the other two starting from 0, joined as the register of the whole chunk: the first's
/// moved on by two thirds, the second's by one, and the third's as it is. The rest goes to
/// [`crc32_arm`].
#[target_feature(enable = "crc,aes")]
fn crc32_arm_3way(mut crc: u32, data: &[u8]) -> u32 {
    let (chunks, rest) = data.as_chunks::<{ 3 * STRIDE }>();
    for chunk in chunks {
        let (a, bc) = chunk.split_at(STRIDE);
        let (b, c) = bc.split_at(STRIDE);
        let (a, b, c) = (
            a.as_chunks::<8>().0,
            b.as_chunks::<8>().0,
            c.as_chunks::<8>().0,
        );
        let (mut ca, mut cb, mut cc) = (crc, 0, 0);
        for ((wa, wb), wc) in a.iter().zip(b).zip(c) {
            ca = __crc32d(ca, u64::from_le_bytes(*wa));
            cb = __crc32d(cb, u64::from_le_bytes(*wb));
            cc = __crc32d(cc, u64::from_le_bytes(*wc));
        }
        crc = shift(ca, SHIFT_2) ^ shift(cb, SHIFT_1) ^ cc;
    }
    crc32_arm(crc, rest)
}

/// The multipliers that fold a 128-bit lane `d` bits forward, for its low and high 64 bits,
/// as in the x86 code.
const fn mults(d: u32) -> [u64; 2] {
    [x_pow_mod_g(d + 31) as u64, x_pow_mod_g(d - 33) as u64]
}

const FOLD_128: [u64; 2] = mults(128);
const FOLD_256: [u64; 2] = mults(256);
const FOLD_384: [u64; 2] = mults(384);
const FOLD_512: [u64; 2] = mults(512);
const FOLD_768: [u64; 2] = mults(768);
const FOLD_1536: [u64; 2] = mults(1536);

#[target_feature(enable = "neon")]
fn load_mults(m: &[u64; 2]) -> poly64x2_t {
    // SAFETY: `m` is 16 readable bytes.
    vreinterpretq_p64_u64(unsafe { vld1q_u64(m.as_ptr()) })
}

#[target_feature(enable = "neon")]
fn load(bytes: &[u8; 16]) -> uint8x16_t {
    // SAFETY: `bytes` is 16 readable bytes; the load needs no alignment.
    unsafe { vld1q_u8(bytes.as_ptr()) }
}

/// `src` moved `d` bits on, by `mults(d)`, and XORed into `dst`.
#[target_feature(enable = "aes,sha3")]
fn fold(src: uint8x16_t, dst: uint8x16_t, mults: poly64x2_t) -> uint8x16_t {
    let src = vreinterpretq_p64_u8(src);
    let low = vmull_p64(vgetq_lane_p64(src, 0), vgetq_lane_p64(mults, 0));
    let high = vmull_high_p64(src, mults);
    veor3q_u8(vreinterpretq_u8_p128(low), vreinterpretq_u8_p128(high), dst)
}

/// The first `16 * N` bytes as `N` accumulators, the register XORed into the first, then
/// each following group of `16 * N` folded on, `mults` moving a lane `16 * N` bytes. Returns
/// the accumulators and the bytes after the last whole group.
#[target_feature(enable = "aes,sha3")]
fn fold_lanes<const N: usize>(
    crc: u32,
    data: &[u8],
    mults: poly64x2_t,
) -> ([uint8x16_t; N], &[u8]) {
    let vecs = data.as_chunks::<16>().0;
    let groups = vecs.len() / N;
    let mut v: [uint8x16_t; N] = std::array::from_fn(|i| load(&vecs[i]));
    v[0] = veorq_u8(
        v[0],
        vreinterpretq_u8_u32(vsetq_lane_u32(crc, vdupq_n_u32(0), 0)),
    );
    for group in vecs[N..groups * N].as_chunks::<N>().0 {
        for i in 0..N {
            v[i] = fold(v[i], load(&group[i]), mults);
        }
    }
    (v, &data[groups * N * 16..])
}

/// The next `N` 16-byte pieces folded into the `N` lanes, `mults` moving a lane `16 * N`
/// bytes, if there are that many; returns the bytes after them.
#[target_feature(enable = "aes,sha3")]
fn fold_pieces<'d, const N: usize>(
    v: &mut [uint8x16_t; N],
    data: &'d [u8],
    mults: poly64x2_t,
) -> &'d [u8] {
    let pieces = data.as_chunks::<16>().0;
    if pieces.len() < N {
        return data;
    }
    for i in 0..N {
        v[i] = fold(v[i], load(&pieces[i]), mults);
    }
    &data[16 * N..]
}

/// Folds with PMULL: twelve lanes for long inputs, four for medium, the CRC instructions
/// for the last bytes and for inputs under 64 bytes.
#[target_feature(enable = "crc,aes,sha3")]
fn crc32_arm_pmull(crc: u32, data: &[u8]) -> u32 {
    if data.len() < 64 {
        return crc32_arm(crc, data);
    }
    let fold_128 = load_mults(&FOLD_128);
    // Each reduction folds a lane onto the one the same distance on, halving the lanes, and
    // while there are as many lanes as pieces left to fill them, folds the pieces in too, so
    // that they fold in parallel rather than one by one.
    let (mut x, rest) = if data.len() < 3 * 192 {
        let (v, rest) = fold_lanes::<4>(crc, data, load_mults(&FOLD_512));
        let m = load_mults(&FOLD_256);
        let mut v = [fold(v[0], v[2], m), fold(v[1], v[3], m)];
        let rest = fold_pieces(&mut v, rest, m);
        v[1] = fold(v[0], v[1], fold_128);
        (v[1], rest)
    } else {
        let (v, rest) = fold_lanes::<12>(crc, data, load_mults(&FOLD_1536));
        let m = load_mults(&FOLD_768);
        let mut v: [_; 6] = std::array::from_fn(|i| fold(v[i], v[i + 6], m));
        let rest = fold_pieces(&mut v, rest, m);
        let m = load_mults(&FOLD_384);
        let mut v: [_; 3] = std::array::from_fn(|i| fold(v[i], v[i + 3], m));
        let rest = fold_pieces(&mut v, rest, m);
        v[1] = fold(v[0], v[1], fold_128);
        (fold(v[1], v[2], fold_128), rest)
    };
    // Fewer pieces left than lanes, folded one by one.
    let (vecs, tail) = rest.as_chunks::<16>();
    for piece in vecs {
        x = fold(x, load(piece), fold_128);
    }
    // The 128 bits reduced to 32 by the CRC instructions, as libdeflate does.
    let x = vreinterpretq_u64_u8(x);
    let crc = __crc32d(0, vgetq_lane_u64(x, 0));
    let crc = __crc32d(crc, vgetq_lane_u64(x, 1));
    crc32_arm(crc, tail)
}

/// The hybrid's vector lanes, rounds per chunk, and words per CRC stream per round: the
/// fastest on an M4 Pro, where the CRC unit then takes about a fifth of the data. More
/// rounds or lanes were slower.
const HYBRID_LANES: usize = 12;
const HYBRID_ROUNDS: usize = 6;
const HYBRID_WORDS: usize = 3;

/// The bytes of a hybrid chunk the vector lanes fold, and each CRC stream's.
const HYBRID_VEC_LEN: usize = 16 * HYBRID_LANES * HYBRID_ROUNDS;
const HYBRID_STREAM_LEN: usize = 8 * HYBRID_WORDS * HYBRID_ROUNDS;
const HYBRID_CHUNK: usize = HYBRID_VEC_LEN + 3 * HYBRID_STREAM_LEN;

/// The multipliers that move a register one, two and three CRC streams on, and a whole
/// chunk on, as [`SHIFT_1`] does for the three-way code.
const HYBRID_SHIFT_1: u32 = x_pow_mod_g(8 * HYBRID_STREAM_LEN as u32 - 33);
const HYBRID_SHIFT_2: u32 = x_pow_mod_g(8 * 2 * HYBRID_STREAM_LEN as u32 - 33);
const HYBRID_SHIFT_3: u32 = x_pow_mod_g(8 * 3 * HYBRID_STREAM_LEN as u32 - 33);
const HYBRID_SHIFT_CHUNK: u32 = x_pow_mod_g(8 * HYBRID_CHUNK as u32 - 33);
const FOLD_HYBRID_LANES: [u64; 2] = mults(128 * HYBRID_LANES as u32);

/// Folds with PMULL and runs the CRC unit at once, on Apple's cores, where they are separate
/// pipes: each chunk is [`HYBRID_LANES`] vector lanes folded over its first
/// [`HYBRID_VEC_LEN`] bytes, then three CRC-instruction streams over three equal parts of
/// the rest, [`HYBRID_WORDS`] words each per round. Every part starts from 0, so chunks don't
/// wait on each other, and joins as in [`crc32_arm_3way`]. What is left after the last whole
/// chunk goes to [`crc32_arm_pmull`].
#[target_feature(enable = "crc,aes,sha3")]
fn crc32_arm_hybrid(mut crc: u32, data: &[u8]) -> u32 {
    let fold_lanes = load_mults(&FOLD_HYBRID_LANES);
    let fold_128 = load_mults(&FOLD_128);
    let (chunks, rest) = data.as_chunks::<HYBRID_CHUNK>();
    for chunk in chunks {
        let (vecs, s) = chunk.split_at(HYBRID_VEC_LEN);
        let (a, s) = s.split_at(HYBRID_STREAM_LEN);
        let (b, c) = s.split_at(HYBRID_STREAM_LEN);
        let groups = vecs.as_chunks::<16>().0.as_chunks::<HYBRID_LANES>().0;
        let (a, b, c) = (
            a.as_chunks::<8>().0.as_chunks::<HYBRID_WORDS>().0,
            b.as_chunks::<8>().0.as_chunks::<HYBRID_WORDS>().0,
            c.as_chunks::<8>().0.as_chunks::<HYBRID_WORDS>().0,
        );
        let mut v: [uint8x16_t; HYBRID_LANES] = std::array::from_fn(|i| load(&groups[0][i]));
        let (mut ca, mut cb, mut cc) = (0, 0, 0);
        for (round, (((group, wa), wb), wc)) in groups.iter().zip(a).zip(b).zip(c).enumerate() {
            if round > 0 {
                for i in 0..HYBRID_LANES {
                    v[i] = fold(v[i], load(&group[i]), fold_lanes);
                }
            }
            for w in 0..HYBRID_WORDS {
                ca = __crc32d(ca, u64::from_le_bytes(wa[w]));
                cb = __crc32d(cb, u64::from_le_bytes(wb[w]));
                cc = __crc32d(cc, u64::from_le_bytes(wc[w]));
            }
        }
        let mut x = v[0];
        for &lane in &v[1..] {
            x = fold(x, lane, fold_128);
        }
        let x = vreinterpretq_u64_u8(x);
        let cv = __crc32d(__crc32d(0, vgetq_lane_u64(x, 0)), vgetq_lane_u64(x, 1));
        crc = shift(crc, HYBRID_SHIFT_CHUNK)
            ^ shift(cv, HYBRID_SHIFT_3)
            ^ shift(ca, HYBRID_SHIFT_2)
            ^ shift(cb, HYBRID_SHIFT_1)
            ^ cc;
    }
    crc32_arm_pmull(crc, rest)
}
