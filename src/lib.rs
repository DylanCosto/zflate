//! A pure-Rust DEFLATE (RFC 1951), zlib (RFC 1950), and gzip (RFC 1952) codec.
//!
//! The API follows libdeflate's: whole buffers in, whole buffers out, no streaming. That is what
//! lets libdeflate skip the bookkeeping a streaming inflater needs, and it is all BGZF and CRAM
//! need.

mod alloc;
mod crc32;
mod deflate;
mod inflate;

// Callers keep a codec per thread and move them between threads: they must stay `Send` and
// `Sync`.
const _: () = {
    const fn send_sync<T: Send + Sync>() {}
    send_sync::<Compressor>();
    send_sync::<Decompressor>();
};

/// Marks the path it is called on as rarely taken, as `std::hint::cold_path` does, which is not
/// stable in Rust 1.91, the oldest that pasta builds with: a call to an empty `#[cold]` function.
#[cold]
#[inline]
pub(crate) fn cold_path() {}

use std::fmt;

pub use crate::crc32::crc32;
#[doc(hidden)]
pub use crate::crc32::implementations as crc32_implementations;
use crate::{
    deflate::Deflater,
    inflate::{Inflated, Tables, inflate, inflate_batch},
};

/// Why a decompression failed, as libdeflate reports it.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum DecompressionError {
    /// The input is not valid for its format, or its checksum or length does not match.
    BadData,
    /// The decompressed data does not fit in the output buffer.
    InsufficientSpace,
}

impl fmt::Display for DecompressionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::BadData => "invalid compressed data",
            Self::InsufficientSpace => "decompressed data does not fit in the output buffer",
        })
    }
}

impl std::error::Error for DecompressionError {}

/// The system could not allocate a compressor's or decompressor's tables: an address-space
/// limit (`ulimit -v`), say. libdeflate's constructors fail the same way.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub struct AllocError;

impl fmt::Display for AllocError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("could not allocate the codec's tables")
    }
}

impl std::error::Error for AllocError {}

/// A decompressor, reusable across calls: it keeps its decode tables, as libdeflate's
/// does.
pub struct Decompressor {
    tables: Box<Tables>,
    /// The second lane's tables, for [`Self::deflate_decompress_batch`], on 64-bit Arm.
    tables_b: Option<Box<Tables>>,
    /// Each stream's result in the batch decoder, kept to reuse its allocation.
    batch_done: Vec<Option<Inflated>>,
}

impl fmt::Debug for Decompressor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Decompressor").finish_non_exhaustive()
    }
}

impl Default for Decompressor {
    fn default() -> Self {
        Self::new()
    }
}

/// gzip header flags (RFC 1952 2.3.1).
const GZIP_FHCRC: u8 = 0x02;
const GZIP_FEXTRA: u8 = 0x04;
const GZIP_FNAME: u8 = 0x08;
const GZIP_FCOMMENT: u8 = 0x10;
const GZIP_FRESERVED: u8 = 0xe0;
/// The CRC32 and ISIZE after a gzip member's data.
const GZIP_FOOTER: usize = 8;

impl Decompressor {
    /// Creates a decompressor.
    ///
    /// # Panics
    ///
    /// If its tables cannot be allocated; [`try_new`](Self::try_new) returns an error instead.
    pub fn new() -> Self {
        Self::try_new().expect("the decompressor's tables can be allocated")
    }

    /// Creates a decompressor, or fails if its tables cannot be allocated.
    ///
    /// # Errors
    ///
    /// [`AllocError`] if the system cannot allocate the tables.
    pub fn try_new() -> Result<Self, AllocError> {
        // The second lane's tables, on 64-bit Arm, where the batch decoder uses them: made here,
        // so that decoding allocates no tables.
        let tables_b = if cfg!(target_arch = "aarch64") {
            Some(Tables::try_new().ok_or(AllocError)?)
        } else {
            None
        };
        Ok(Self {
            tables: Tables::try_new().ok_or(AllocError)?,
            tables_b,
            batch_done: Vec::new(),
        })
    }

    /// Inflates raw DEFLATE data into `out`, returning the length written. Input after the
    /// final block is ignored.
    ///
    /// # Errors
    ///
    /// [`DecompressionError::BadData`] if the data is not valid DEFLATE, and
    /// [`DecompressionError::InsufficientSpace`] if it inflates to more than `out` holds.
    pub fn deflate_decompress(
        &mut self,
        data: &[u8],
        out: &mut [u8],
    ) -> Result<usize, DecompressionError> {
        Ok(inflate(&mut self.tables, data, out)?.written)
    }

    /// Inflates each of `data`, independent raw DEFLATE streams, into the same-numbered
    /// buffer of `out`, and writes each stream's result to the same-numbered element of
    /// `results`: the BGZF blocks of a file, for example. Each result is what
    /// [`deflate_decompress`](Self::deflate_decompress) would return for that stream.
    ///
    /// On 64-bit Arm it decodes two streams at a time, so that a wide core overlaps their
    /// work; elsewhere it decodes them one by one. A stream that fails in a pair is decoded
    /// again alone, for exactly the single decoder's error; as with that decoder, a failed
    /// stream's buffer then holds unspecified bytes. The first batch allocates a result slot
    /// per stream, reused by later batches of the same size or smaller.
    ///
    /// # Errors
    ///
    /// Each stream's error is in its element of `results`, and does not stop the others:
    /// [`DecompressionError::BadData`] if it is not valid DEFLATE, and
    /// [`DecompressionError::InsufficientSpace`] if it inflates to more than its buffer
    /// holds.
    ///
    /// # Panics
    ///
    /// If `data`, `out` and `results` differ in length.
    pub fn deflate_decompress_batch(
        &mut self,
        data: &[&[u8]],
        out: &mut [&mut [u8]],
        results: &mut [Result<usize, DecompressionError>],
    ) {
        assert!(
            data.len() == out.len() && data.len() == results.len(),
            "one output buffer and one result per stream: {} streams, {} buffers, {} results",
            data.len(),
            out.len(),
            results.len()
        );
        let mut done = std::mem::take(&mut self.batch_done);
        done.clear();
        done.resize(data.len(), None);
        if let Some(tables_b) = &mut self.tables_b {
            inflate_batch([&mut self.tables, tables_b], data, &mut *out, &mut done);
        }
        for (((r, d), o), result) in done.iter().zip(data).zip(out.iter_mut()).zip(results) {
            *result = match r {
                Some(i) => Ok(i.written),
                // Failed in a lane, or not tried: alone, for the single decoder's error.
                None => self.deflate_decompress(d, o),
            };
        }
        self.batch_done = done;
    }

    /// Inflates zlib data into `out`, checking its Adler-32, and returns the length written.
    ///
    /// # Errors
    ///
    /// As [`deflate_decompress`](Self::deflate_decompress), and
    /// [`DecompressionError::BadData`] if the header or the Adler-32 is wrong.
    pub fn zlib_decompress(
        &mut self,
        data: &[u8],
        out: &mut [u8],
    ) -> Result<usize, DecompressionError> {
        use DecompressionError::BadData;
        // CMF and FLG, then the data, then the Adler-32.
        let [cmf, flg, ..] = *data else {
            return Err(BadData);
        };
        let body = data.get(2..data.len().wrapping_sub(4)).ok_or(BadData)?;
        let header = u16::from_be_bytes([cmf, flg]);
        // FCHECK, CM 8 (DEFLATE), CINFO at most 7 (32 KiB window), and no FDICT.
        if header % 31 != 0 || cmf & 0xf != 8 || cmf >> 4 > 7 || flg & 0x20 != 0 {
            return Err(BadData);
        }
        let Inflated { consumed, written } = inflate(&mut self.tables, body, out)?;
        let trailer = &data[2 + consumed..][..4];
        if adler32(1, &out[..written]).to_be_bytes() != *trailer {
            return Err(BadData);
        }
        Ok(written)
    }

    /// Inflates one gzip member into `out`, checking its CRC32 and length, and returns the
    /// length written. Input after the member is ignored.
    ///
    /// # Errors
    ///
    /// As [`deflate_decompress`](Self::deflate_decompress), and
    /// [`DecompressionError::BadData`] if the header, the CRC32, or the length is wrong.
    pub fn gzip_decompress(
        &mut self,
        data: &[u8],
        out: &mut [u8],
    ) -> Result<usize, DecompressionError> {
        use DecompressionError::BadData;
        // ID1, ID2, CM, FLG, MTIME, XFL, OS: 10 bytes, and at least the footer after them.
        if data.len() < 10 + GZIP_FOOTER || data[..3] != [0x1f, 0x8b, 8] {
            return Err(BadData);
        }
        let flg = data[3];
        if flg & GZIP_FRESERVED != 0 {
            return Err(BadData);
        }
        let room = |pos: usize| pos + GZIP_FOOTER <= data.len();
        let mut pos = 10;
        if flg & GZIP_FEXTRA != 0 {
            let xlen = usize::from(u16::from_le_bytes([data[10], data[11]]));
            pos += 2 + xlen;
            if !room(pos) {
                return Err(BadData);
            }
        }
        // The file name and comment are zero-terminated; like libdeflate, stop at the end
        // of the input if there is no terminator.
        for field in [GZIP_FNAME, GZIP_FCOMMENT] {
            if flg & field != 0 {
                while data[pos] != 0 && pos + 1 < data.len() {
                    pos += 1;
                }
                pos += 1;
                if !room(pos) {
                    return Err(BadData);
                }
            }
        }
        if flg & GZIP_FHCRC != 0 {
            pos += 2;
            if !room(pos) {
                return Err(BadData);
            }
        }
        let body = &data[pos..data.len() - GZIP_FOOTER];
        let Inflated { consumed, written } = inflate(&mut self.tables, body, out)?;
        let footer = &data[pos + consumed..][..GZIP_FOOTER];
        if crc32(0, &out[..written]).to_le_bytes() != footer[..4]
            || (written as u32).to_le_bytes() != footer[4..]
        {
            return Err(BadData);
        }
        Ok(written)
    }
}

/// A compression level, as libdeflate's: 0 stores the data uncompressed, 1 is the fastest,
/// and 12 compresses the most.
///
/// Each level compresses as libdeflate's does, to within about 0.05% in size, except level
/// 5, which searches less deeply and takes no 3-byte matches: on BAM, BCF, VCF and FASTA
/// files it is 21-77% faster than libdeflate's level 5, and its output is between 1% smaller
/// and 1.5% larger. Other data can grow more (CRAM's gzip blocks of unaligned HiFi reads:
/// 8.8%).
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct CompressionLevel(u8);

impl CompressionLevel {
    /// No compression: stored blocks.
    pub const NONE: Self = Self(0);
    /// The fastest compression.
    pub const FASTEST: Self = Self(1);
    /// libdeflate's (and zlib's) default.
    pub const DEFAULT: Self = Self(6);
    /// The most compression.
    pub const BEST: Self = Self(12);

    /// The level `n`, if it is 0 to 12.
    pub const fn new(n: u8) -> Option<Self> {
        if n <= 12 { Some(Self(n)) } else { None }
    }

    /// The level as a number.
    pub const fn get(self) -> u8 {
        self.0
    }
}

impl Default for CompressionLevel {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Why a compression failed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum CompressionError {
    /// The compressed data does not fit in the output buffer. A buffer of the
    /// `*_compress_bound` length always fits.
    InsufficientSpace,
}

impl fmt::Display for CompressionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InsufficientSpace => "compressed data does not fit in the output buffer",
        })
    }
}

impl std::error::Error for CompressionError {}

/// A compressor at one level, reusable across calls: it keeps its match finder's tables,
/// as libdeflate's does. The output need not match libdeflate's byte for byte.
pub struct Compressor {
    level: CompressionLevel,
    deflater: Box<Deflater>,
}

impl fmt::Debug for Compressor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Compressor")
            .field("level", &self.level)
            .finish_non_exhaustive()
    }
}

/// zlib's and gzip's headers and trailers.
const ZLIB_OVERHEAD: usize = 2 + 4;
const GZIP_HEADER: usize = 10;
const GZIP_OVERHEAD: usize = GZIP_HEADER + GZIP_FOOTER;

impl Compressor {
    /// Creates a compressor at `level`.
    ///
    /// # Panics
    ///
    /// If its tables cannot be allocated; [`try_new`](Self::try_new) returns an error instead.
    pub fn new(level: CompressionLevel) -> Self {
        Self::try_new(level).expect("the compressor's tables can be allocated")
    }

    /// Creates a compressor at `level`, or fails if its tables cannot be allocated. Levels
    /// 10 to 12 also grow buffers as they compress, with ordinary allocations.
    ///
    /// # Errors
    ///
    /// [`AllocError`] if the system cannot allocate the tables.
    pub fn try_new(level: CompressionLevel) -> Result<Self, AllocError> {
        let deflater = Deflater::try_new(level.get()).and_then(alloc::try_box);
        Ok(Self {
            level,
            deflater: deflater.ok_or(AllocError)?,
        })
    }

    /// The compressor's level.
    pub fn level(&self) -> CompressionLevel {
        self.level
    }

    /// Compresses `data` into `out` as raw DEFLATE, returning the length written.
    ///
    /// # Errors
    ///
    /// [`CompressionError::InsufficientSpace`] if the result does not fit in `out`, which
    /// it always does with [`deflate_compress_bound`] bytes.
    pub fn deflate_compress(
        &mut self,
        data: &[u8],
        out: &mut [u8],
    ) -> Result<usize, CompressionError> {
        self.deflater
            .compress(data, out)
            .ok_or(CompressionError::InsufficientSpace)
    }

    /// Compresses `data` into `out` as a zlib stream, returning the length written.
    ///
    /// # Errors
    ///
    /// As [`deflate_compress`](Self::deflate_compress), with [`zlib_compress_bound`].
    pub fn zlib_compress(
        &mut self,
        data: &[u8],
        out: &mut [u8],
    ) -> Result<usize, CompressionError> {
        if out.len() <= ZLIB_OVERHEAD {
            return Err(CompressionError::InsufficientSpace);
        }
        // CMF: DEFLATE with a 32 KiB window; FLG: the level's hint, and the check bits that
        // make the header a multiple of 31.
        let hint: u16 = match self.level.get() {
            0..2 => 0,
            2..6 => 1,
            6 => 2,
            _ => 3,
        };
        let mut header = 0x7800 | hint << 6;
        header |= 31 - header % 31;
        out[..2].copy_from_slice(&header.to_be_bytes());
        let end = out.len() - 4;
        let len = self.deflate_compress(data, &mut out[2..end])?;
        out[2 + len..][..4].copy_from_slice(&adler32(1, data).to_be_bytes());
        Ok(2 + len + 4)
    }

    /// Compresses `data` into `out` as a gzip member, returning the length written. The
    /// header has no name, time or extra field.
    ///
    /// # Errors
    ///
    /// As [`deflate_compress`](Self::deflate_compress), with [`gzip_compress_bound`].
    pub fn gzip_compress(
        &mut self,
        data: &[u8],
        out: &mut [u8],
    ) -> Result<usize, CompressionError> {
        if out.len() <= GZIP_OVERHEAD {
            return Err(CompressionError::InsufficientSpace);
        }
        // XFL tells the fastest and the slowest levels apart; the OS is unknown.
        let xfl = match self.level.get() {
            0..2 => 4,
            2..8 => 0,
            _ => 2,
        };
        out[..GZIP_HEADER].copy_from_slice(&[0x1f, 0x8b, 8, 0, 0, 0, 0, 0, xfl, 0xff]);
        let end = out.len() - GZIP_FOOTER;
        let len = self.deflate_compress(data, &mut out[GZIP_HEADER..end])?;
        let footer = &mut out[GZIP_HEADER + len..][..GZIP_FOOTER];
        footer[..4].copy_from_slice(&crc32(0, data).to_le_bytes());
        footer[4..].copy_from_slice(&(data.len() as u32).to_le_bytes());
        Ok(GZIP_HEADER + len + GZIP_FOOTER)
    }
}

/// The most raw DEFLATE data `len` bytes can compress to, at any level.
pub fn deflate_compress_bound(len: usize) -> usize {
    deflate::compress_bound(len)
}

/// The most a zlib stream of `len` bytes can compress to, at any level.
pub fn zlib_compress_bound(len: usize) -> usize {
    ZLIB_OVERHEAD + deflate_compress_bound(len)
}

/// The most a gzip member of `len` bytes can compress to, at any level.
pub fn gzip_compress_bound(len: usize) -> usize {
    GZIP_OVERHEAD + deflate_compress_bound(len)
}

/// Updates an Adler-32 with `data`, starting from 1, as libdeflate's `libdeflate_adler32` does.
pub fn adler32(adler: u32, data: &[u8]) -> u32 {
    let mut hash = simd_adler32::Adler32::from_checksum(adler);
    hash.write(data);
    hash.finish()
}
