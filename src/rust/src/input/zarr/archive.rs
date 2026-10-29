//! Minimal zarr v2 reader for a zip archive of STORED, blosc-compressed
//! chunks: the layout of an Atera v1 `transcripts.zarr.zip`.

use std::fs::File;
use std::io::Read;

use serde_json::Value;
use zip::result::ZipError;
use zip::ZipArchive;

pub struct Archive {
    zip: ZipArchive<File>,
}

pub struct ArrayMeta {
    pub shape: Vec<usize>,
    pub chunks: Vec<usize>,
    pub dtype: String,
    pub order_f: bool,
}

impl ArrayMeta {
    pub fn ncol(&self) -> usize {
        if self.shape.len() > 1 { self.shape[1] } else { 1 }
    }
}

impl Archive {
    pub fn open(path: &str) -> Result<Self, String> {
        let f = File::open(path).map_err(|e| format!("open {path}: {e}"))?;
        let zip = ZipArchive::new(f).map_err(|e| format!("zip {path}: {e}"))?;
        Ok(Archive { zip })
    }

    pub fn read(&mut self, name: &str) -> Result<Option<Vec<u8>>, String> {
        match self.zip.by_name(name) {
            Ok(mut f) => {
                let mut v = Vec::with_capacity(f.size() as usize);
                f.read_to_end(&mut v).map_err(|e| format!("read {name}: {e}"))?;
                Ok(Some(v))
            }
            Err(ZipError::FileNotFound) => Ok(None),
            Err(e) => Err(format!("zip entry {name}: {e}")),
        }
    }

    pub fn meta(&mut self, prefix: &str) -> Result<ArrayMeta, String> {
        let name = format!("{prefix}/.zarray");
        let bytes = self.read(&name)?.ok_or(format!("missing {name}"))?;
        let v: Value = serde_json::from_slice(&bytes).map_err(|e| format!("{name}: {e}"))?;
        let ints = |k: &str| -> Result<Vec<usize>, String> {
            v[k].as_array()
                .ok_or(format!("{name}: no {k}"))?
                .iter()
                .map(|x| x.as_u64().map(|u| u as usize).ok_or(format!("{name}: bad {k}")))
                .collect()
        };
        Ok(ArrayMeta {
            shape: ints("shape")?,
            chunks: ints("chunks")?,
            dtype: v["dtype"].as_str().unwrap_or("").to_string(),
            order_f: v["order"].as_str() == Some("F"),
        })
    }

    /// Decompressed bytes of row-chunk `ci` (the `<ci>.0` key).
    fn chunk(&mut self, prefix: &str, ci: usize) -> Result<Vec<u8>, String> {
        let name = format!("{prefix}/{ci}.0");
        // fill_value is null in these archives, so a missing chunk would be
        // NA in the R reader; refuse rather than guess
        let raw = self.read(&name)?.ok_or(format!("missing chunk {name}"))?;
        blosc_decompress(&raw).map_err(|e| format!("{name}: {e}"))
    }

    /// Rows `[r0, r1)` of every column of a 2D array, decoded by `conv`
    /// from `size`-byte little-endian elements.
    pub fn read_cols<T: Copy>(
        &mut self,
        prefix: &str,
        meta: &ArrayMeta,
        r0: usize,
        r1: usize,
        size: usize,
        conv: impl Fn(&[u8]) -> T,
    ) -> Result<Vec<Vec<T>>, String> {
        let nc = meta.ncol();
        let cn = meta.chunks[0];
        let n = r1 - r0;
        let mut out: Vec<Vec<T>> = (0..nc).map(|_| Vec::with_capacity(n)).collect();
        if n == 0 {
            return Ok(out);
        }
        for ci in (r0 / cn)..=((r1 - 1) / cn) {
            let bytes = self.chunk(prefix, ci)?;
            if bytes.len() < cn * nc * size {
                return Err(format!("{prefix}/{ci}.0: short chunk"));
            }
            let cs = ci * cn;
            let a = r0.max(cs) - cs;
            let b = r1.min(cs + cn) - cs;
            for (j, col) in out.iter_mut().enumerate() {
                for i in a..b {
                    let k = if meta.order_f { j * cn + i } else { i * nc + j };
                    col.push(conv(&bytes[k * size..(k + 1) * size]));
                }
            }
        }
        Ok(out)
    }
}

fn blosc_decompress(src: &[u8]) -> Result<Vec<u8>, String> {
    let mut nbytes = 0usize;
    let mut cbytes = 0usize;
    let mut blocksize = 0usize;
    unsafe {
        blosc_src::blosc_cbuffer_sizes(
            src.as_ptr() as *const _, &mut nbytes, &mut cbytes, &mut blocksize,
        );
    }
    if cbytes > src.len() {
        return Err("truncated blosc buffer".into());
    }
    let mut out = vec![0u8; nbytes];
    let got = unsafe {
        blosc_src::blosc_decompress_ctx(
            src.as_ptr() as *const _, out.as_mut_ptr() as *mut _, nbytes, 1,
        )
    };
    if got < 0 || got as usize != nbytes {
        return Err(format!("blosc decompress failed ({got})"));
    }
    Ok(out)
}

pub fn f32_le(b: &[u8]) -> f32 { f32::from_le_bytes([b[0], b[1], b[2], b[3]]) }
pub fn u32_le(b: &[u8]) -> u32 { u32::from_le_bytes([b[0], b[1], b[2], b[3]]) }
pub fn f16_le(b: &[u8]) -> f32 { half::f16::from_le_bytes([b[0], b[1]]).to_f32() }
pub fn u8_le(b: &[u8]) -> u8 { b[0] }
