//! The zlib Adler-32 (RFC 1950): two sums modulo 65521, `s1` of the bytes and `s2` of the
//! running `s1`s, in one 32-bit value, `s2 << 16 | s1`.
//!
//! The fast paths are ports of libdeflate's `x86/adler32_template.h` and `arm/adler32_impl.h`
//! (MIT); [`scalar`] is the fallback and handles what the vector code leaves at the end of
//! each chunk. The best implementation for the CPU is chosen once, at the first call.

#[cfg(target_arch = "aarch64")]
mod arm;
mod scalar;
#[cfg(target_arch = "x86_64")]
mod x86;

use std::sync::OnceLock;

/// An Adler-32 implementation. The checksum it starts from has both halves below [`DIVISOR`].
pub type Impl = fn(u32, &[u8]) -> u32;

/// The modulus of both sums, the largest prime below 2^16.
const DIVISOR: u32 = 65_521;

/// The most bytes that can be added before `s2` could overflow 32 bits, when both sums start
/// at `DIVISOR - 1` and every byte is `0xff` (libdeflate's `MAX_CHUNK_LEN`).
const MAX_CHUNK_LEN: usize = 5552;

/// Updates an Adler-32 with `data`, starting from 1, as libdeflate's `libdeflate_adler32` does.
///
/// `adler32(adler32(1, a), b)` is the Adler-32 of `a` followed by `b`. Halves at or above
/// 65521, which no data checksums to, are reduced modulo 65521 first.
#[must_use]
pub fn adler32(adler: u32, data: &[u8]) -> u32 {
    if data.is_empty() {
        return adler;
    }
    let adler = (((adler >> 16) % DIVISOR) << 16) | ((adler & 0xffff) % DIVISOR);
    best()(adler, data)
}

fn best() -> Impl {
    static BEST: OnceLock<Impl> = OnceLock::new();
    *BEST.get_or_init(|| available().first().map_or(scalar::generic, |&(_, f)| f))
}

/// The implementations this CPU supports, fastest first.
fn available() -> Vec<(&'static str, Impl)> {
    let mut impls = Vec::new();
    #[cfg(target_arch = "x86_64")]
    if !cfg!(feature = "force-portable") {
        impls.extend(x86::available());
    }
    // Not under Miri, which cannot run the inline assembly of the dot-product version.
    #[cfg(target_arch = "aarch64")]
    if !cfg!(feature = "force-portable") && !cfg!(miri) {
        impls.extend(arm::available());
    }
    impls.push(("generic", scalar::generic as Impl));
    impls
}

/// Each Adler-32 implementation this CPU supports, by name, fastest first, for tests and
/// benchmarks. Not part of the stable API.
#[doc(hidden)]
pub fn implementations() -> Vec<(&'static str, Impl)> {
    available()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The definition, a byte at a time.
    fn reference(adler: u32, data: &[u8]) -> u32 {
        let (mut s1, mut s2) = ((adler & 0xffff) % DIVISOR, (adler >> 16) % DIVISOR);
        for &b in data {
            s1 = (s1 + u32::from(b)) % DIVISOR;
            s2 = (s2 + s1) % DIVISOR;
        }
        (s2 << 16) | s1
    }

    /// Pseudorandom bytes.
    fn noise(len: usize) -> Vec<u8> {
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (state >> 24) as u8
            })
            .collect()
    }

    /// Starting checksums: the initial one, zero, and both sums at their largest.
    const STARTS: [u32; 3] = [1, 0, ((DIVISOR - 1) << 16) | (DIVISOR - 1)];

    #[test]
    fn check_value() {
        for (name, f) in available() {
            assert_eq!(f(1, b"Wikipedia"), 0x11e6_0398, "{name}");
        }
    }

    /// Every length up to past four 512-bit vectors twice over, and the lengths around each
    /// chunk limit, with all-`0xff` bytes (the largest sums) and pseudorandom ones.
    #[test]
    fn matches_the_definition() {
        let mut lens: Vec<usize> = (0..=1100).collect();
        for edge in [2048, 4096, 5504, 5552, 8192, 11_008, 11_104, 65_536] {
            lens.extend(edge - 3..=edge + 3);
        }
        lens.push(1 << 20);
        let max = *lens.iter().max().expect("lengths");

        for data in [vec![0xff; max], noise(max)] {
            for (name, f) in available() {
                for &len in &lens {
                    for start in STARTS {
                        let s = &data[..len];
                        assert_eq!(
                            f(start, s),
                            reference(start, s),
                            "{name}, {len} bytes from {start:#x}"
                        );
                    }
                }
            }
        }
    }

    /// Long inputs that do not start on a vector boundary, which the vector code aligns first.
    #[test]
    fn misaligned_long_inputs() {
        let data = noise(70_000 + 64);
        let ff = vec![0xff; 70_000 + 64];
        for (name, f) in available() {
            for offset in 0..64 {
                for d in [&data, &ff] {
                    let s = &d[offset..offset + 70_000];
                    assert_eq!(f(1, s), reference(1, s), "{name}, offset {offset}");
                }
            }
        }
    }

    #[test]
    fn unreduced_starts_and_chaining() {
        let data = noise(100_000);
        assert_eq!(adler32(0xffff_ffff, &data), reference(0xffff_ffff, &data));
        assert_eq!(adler32(0xffff_ffff, &[]), 0xffff_ffff);
        for split in [0, 1, 63, 64, 5552, 50_000, 100_000] {
            let (a, b) = data.split_at(split);
            assert_eq!(
                adler32(adler32(1, a), b),
                adler32(1, &data),
                "split at {split}"
            );
        }
    }
}
