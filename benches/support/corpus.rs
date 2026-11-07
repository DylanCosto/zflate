//! Real BGZF data for tests and benchmarks: the blocks of BAM, BCF, and bgzipped files, which
//! is what BGZF readers and writers hand a DEFLATE codec.
//!
//! `ZFLATE_CORPUS` lists the files, separated by `:`. Without it there are none: the parts of
//! the tests that need real data are skipped (the differential tests say so), and the
//! benchmarks have nothing to run on (the examples take a file instead). Files that do not
//! exist are skipped. A file larger than `ZFLATE_WINDOW_MB` (default 48) is sampled: that
//! many megabytes of whole blocks from its middle, which is more typical of it than its
//! start. `path@0.3` centres the window at 30% of the file instead, so one large file can
//! give separate windows for tuning and checking.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// One BGZF block: its raw DEFLATE data and the CRC32 and length of what it inflates to.
pub struct Block {
    pub deflate: Vec<u8>,
    pub crc32: u32,
    pub len: usize,
}

/// A corpus file split into its BGZF blocks.
pub struct File {
    pub name: String,
    pub blocks: Vec<Block>,
}

impl File {
    /// The total length of the blocks' inflated data.
    pub fn inflated_len(&self) -> usize {
        self.blocks.iter().map(|b| b.len).sum()
    }

    /// The total length of the blocks' DEFLATE data.
    pub fn deflated_len(&self) -> usize {
        self.blocks.iter().map(|b| b.deflate.len()).sum()
    }
}

/// The corpus files that exist.
pub fn files() -> Vec<File> {
    paths().iter().filter_map(|p| read(p)).collect()
}

fn paths() -> Vec<PathBuf> {
    std::env::var("ZFLATE_CORPUS")
        .unwrap_or_default()
        .split(':')
        .filter(|s| !s.is_empty())
        .map(PathBuf::from)
        .collect()
}

fn read(spec: &Path) -> Option<File> {
    let spec = spec.to_str()?;
    let (path, at) = match spec.rsplit_once('@') {
        Some((p, f)) => (Path::new(p), f.parse().ok()?),
        None => (Path::new(spec), 0.5),
    };
    let data = window(path, at)?;
    let blocks = split(&data).unwrap_or_else(|| panic!("{} is not BGZF", path.display()));
    let name = Path::new(spec).file_name()?.to_string_lossy().into_owned();
    Some(File { name, blocks })
}

/// The whole file, or for a large one, whole blocks from a window centred at fraction `at`
/// of it.
fn window(path: &Path, at: f64) -> Option<Vec<u8>> {
    use std::io::{Read, Seek, SeekFrom};

    let limit: u64 = std::env::var("ZFLATE_WINDOW_MB")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(48)
        << 20;
    let mut f = std::fs::File::open(path).ok()?;
    let size = f.metadata().ok()?.len();
    if size <= limit {
        return std::fs::read(path).ok();
    }
    let centre = (size as f64 * at) as u64;
    f.seek(SeekFrom::Start(
        centre.clamp(limit / 2, size - limit / 2) - limit / 2,
    ))
    .ok()?;
    let mut buf = vec![0; limit as usize];
    f.read_exact(&mut buf).ok()?;
    // The first block start whose next block also starts where its size says.
    let start = (0..buf.len() - 65_536).find(|&i| {
        let at = |i: usize| block_size(&buf[i..]);
        at(i).is_some_and(|n| at(i + n).is_some())
    })?;
    // Up to the last whole block.
    let mut end = start;
    while let Some(n) = buf.get(end..).and_then(block_size) {
        if end + n > buf.len() {
            break;
        }
        end += n;
    }
    Some(buf[start..end].to_vec())
}

/// The size of the BGZF block at the start of `data`, if one starts there.
fn block_size(data: &[u8]) -> Option<usize> {
    if data.get(..4)? != [0x1f, 0x8b, 8, 4] {
        return None;
    }
    let xlen = usize::from(u16::from_le_bytes(data.get(10..12)?.try_into().ok()?));
    bc_subfield(data.get(12..12 + xlen)?).map(|b| b + 1)
}

/// Splits BGZF data into its blocks, or `None` if it is not BGZF.
fn split(mut data: &[u8]) -> Option<Vec<Block>> {
    let mut blocks = Vec::new();
    while !data.is_empty() {
        // ID1 ID2 CM FLG(FEXTRA), MTIME, XFL, OS, XLEN, then the BC subfield with BSIZE.
        if data.get(..4)? != [0x1f, 0x8b, 8, 4] {
            return None;
        }
        let xlen = usize::from(u16::from_le_bytes(data.get(10..12)?.try_into().ok()?));
        let extra = data.get(12..12 + xlen)?;
        let bsize = bc_subfield(extra)?;
        let block = data.get(..bsize + 1)?;
        let trailer = block.get(block.len() - 8..)?;
        blocks.push(Block {
            deflate: block.get(12 + xlen..block.len() - 8)?.to_vec(),
            crc32: u32::from_le_bytes(trailer[..4].try_into().ok()?),
            len: u32::from_le_bytes(trailer[4..].try_into().ok()?) as usize,
        });
        data = &data[block.len()..];
    }
    Some(blocks)
}

fn bc_subfield(mut extra: &[u8]) -> Option<usize> {
    while extra.len() >= 4 {
        let len = usize::from(u16::from_le_bytes([extra[2], extra[3]]));
        if extra[..2] == *b"BC" && len == 2 {
            return Some(usize::from(u16::from_le_bytes(
                extra.get(4..6)?.try_into().ok()?,
            )));
        }
        extra = extra.get(4 + len..)?;
    }
    None
}
