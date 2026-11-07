//! Table-driven CRC32: the fallback without carry-less multiplication, and the tail of
//! inputs too short to fold.

use super::POLY;

/// `TABLES[k][b]` is the CRC update for byte `b` followed by `k` zero bytes, so eight bytes
/// can be folded in with eight independent lookups (slice-by-8).
static TABLES: [[u32; 256]; 8] = tables();

const fn tables() -> [[u32; 256]; 8] {
    let mut t = [[0; 256]; 8];
    let mut b = 0;
    while b < 256 {
        let mut crc = b as u32;
        let mut bit = 0;
        while bit < 8 {
            crc = (crc >> 1) ^ if crc & 1 != 0 { POLY } else { 0 };
            bit += 1;
        }
        t[0][b] = crc;
        b += 1;
    }
    let mut k = 1;
    while k < 8 {
        let mut b = 0;
        while b < 256 {
            let prev = t[k - 1][b];
            t[k][b] = (prev >> 8) ^ t[0][(prev & 0xff) as usize];
            b += 1;
        }
        k += 1;
    }
    t
}

/// One byte at a time; the shortest code, for a few bytes.
pub(super) fn slice1(mut crc: u32, data: &[u8]) -> u32 {
    for &b in data {
        crc = (crc >> 8) ^ TABLES[0][usize::from(crc as u8 ^ b)];
    }
    crc
}

/// Eight bytes at a time.
pub(super) fn slice8(mut crc: u32, data: &[u8]) -> u32 {
    let (chunks, rest) = data.as_chunks::<8>();
    for chunk in chunks {
        let v = u64::from_le_bytes(*chunk) ^ u64::from(crc);
        let byte = |i: u32| usize::from((v >> (8 * i)) as u8);
        crc = TABLES[7][byte(0)]
            ^ TABLES[6][byte(1)]
            ^ TABLES[5][byte(2)]
            ^ TABLES[4][byte(3)]
            ^ TABLES[3][byte(4)]
            ^ TABLES[2][byte(5)]
            ^ TABLES[1][byte(6)]
            ^ TABLES[0][byte(7)];
    }
    slice1(crc, rest)
}
