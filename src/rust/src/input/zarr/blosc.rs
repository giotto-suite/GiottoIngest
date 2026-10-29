//! Blosc1 decompression in Rust, for the codecs Atera v1 archives use (zstd,
//! byte shuffle). Follows c-blosc 1.x `blosc_run_decompression_with_context`
//! and `blosc_d`. c-blosc itself is not linked: on Windows it compiles its
//! own pthread implementation, which clashes with MinGW's winpthread.

const MIN_HEADER: usize = 16;
const VERSION_FORMAT: u8 = 2;
const DOSHUFFLE: u8 = 0x1;
const MEMCPYED: u8 = 0x2;
const DOBITSHUFFLE: u8 = 0x4;
const DONT_SPLIT: u8 = 0x10;
const ZSTD_FORMAT: u8 = 4;
const MAX_SPLITS: usize = 16;
const MIN_BUFFERSIZE: usize = 128;

fn u32_at(b: &[u8], at: usize) -> Result<usize, String> {
    b.get(at..at + 4)
        .map(|s| u32::from_le_bytes([s[0], s[1], s[2], s[3]]) as usize)
        .ok_or_else(|| "truncated blosc buffer".to_string())
}

/// Decompresses one Blosc1 frame.
pub fn decompress(src: &[u8]) -> Result<Vec<u8>, String> {
    if src.len() < MIN_HEADER {
        return Err("truncated blosc buffer".into());
    }
    let (version, flags, typesize) = (src[0], src[2], src[3] as usize);
    let nbytes = u32_at(src, 4)?;
    let blocksize = u32_at(src, 8)?;
    let cbytes = u32_at(src, 12)?;
    if nbytes == 0 {
        return Ok(Vec::new());
    }
    if version != VERSION_FORMAT || flags & 0x08 != 0 {
        return Err(format!("unsupported: blosc format version {version}, flags {flags:#x}"));
    }
    if cbytes > src.len() || blocksize == 0 || typesize == 0 {
        return Err("corrupt blosc header".into());
    }
    let mut out = vec![0u8; nbytes];
    if flags & MEMCPYED != 0 {
        let body = src.get(MIN_HEADER..MIN_HEADER + nbytes).ok_or("truncated blosc buffer")?;
        out.copy_from_slice(body);
        return Ok(out);
    }
    let codec = (flags & 0xe0) >> 5;
    if codec != ZSTD_FORMAT {
        return Err(format!("unsupported: blosc codec {codec} (only zstd)"));
    }
    if flags & DOBITSHUFFLE != 0 && blocksize >= typesize {
        return Err("unsupported: blosc bitshuffle".into());
    }
    let dont_split = flags & DONT_SPLIT != 0;
    let nblocks = nbytes.div_ceil(blocksize);
    let leftover = nbytes % blocksize;
    let mut tmp = vec![0u8; blocksize];
    for j in 0..nblocks {
        let last = j == nblocks - 1 && leftover > 0;
        let bsize = if last { leftover } else { blocksize };
        let doshuffle = flags & DOSHUFFLE != 0 && typesize > 1;
        let dest = &mut out[j * blocksize..j * blocksize + bsize];
        let target: &mut [u8] = if doshuffle { &mut tmp[..bsize] } else { dest };
        let nsplits = if !dont_split && typesize <= MAX_SPLITS && bsize / typesize >= MIN_BUFFERSIZE
            && !last { typesize } else { 1 };
        let neblock = bsize / nsplits;
        let mut at = u32_at(src, MIN_HEADER + j * 4)?;
        for s in 0..nsplits {
            let n = u32_at(src, at)?;
            at += 4;
            let part = src.get(at..at + n).ok_or("truncated blosc block")?;
            let into = &mut target[s * neblock..(s + 1) * neblock];
            if n == neblock {
                into.copy_from_slice(part);
            } else {
                let got = zstd::bulk::decompress_to_buffer(part, into)
                    .map_err(|e| format!("blosc zstd block: {e}"))?;
                if got != neblock {
                    return Err(format!("blosc block: {got} bytes where {neblock} expected"));
                }
            }
            at += n;
        }
        if doshuffle {
            unshuffle(typesize, &tmp[..bsize], &mut out[j * blocksize..j * blocksize + bsize]);
        }
    }
    Ok(out)
}

/// Byte unshuffle of one block; trailing bytes past whole elements are
/// copied as they are.
fn unshuffle(ts: usize, src: &[u8], dest: &mut [u8]) {
    let ne = src.len() / ts;
    for i in 0..ne {
        for k in 0..ts {
            dest[i * ts + k] = src[k * ne + i];
        }
    }
    let rem = ne * ts;
    dest[rem..].copy_from_slice(&src[rem..]);
}

#[cfg(test)]
mod tests {
    use std::io::Read;

    /// Every chunk of the zip archives in `GI_BLOSC_PARITY` (colon separated)
    /// decodes to the same bytes as c-blosc. Without the variable, the tiny
    /// Atera dataset under GI_DATA (default ~/Documents/Datasets), if present.
    #[test]
    fn matches_c_blosc() {
        let paths: Vec<String> = match std::env::var("GI_BLOSC_PARITY") {
            Ok(v) => v.split(':').map(String::from).collect(),
            Err(_) => {
                let d = std::env::var("GI_DATA").unwrap_or_else(|_| {
                    format!("{}/Documents/Datasets", std::env::var("HOME").unwrap_or_default())
                });
                vec![format!("{d}/Atera/tiny_atera_dataset_outs/transcripts.zarr.zip")]
            }
        };
        for p in paths.iter().filter(|p| std::path::Path::new(p).exists()) {
            let mut z = zip::ZipArchive::new(std::fs::File::open(p).unwrap()).unwrap();
            let names: Vec<String> = z.file_names().map(String::from).collect();
            let mut n = 0;
            for name in names.iter().filter(|n| {
                let last = n.rsplit('/').next().unwrap_or("");
                !last.starts_with('.') && !n.ends_with('/')
            }) {
                let mut raw = Vec::new();
                z.by_name(name).unwrap().read_to_end(&mut raw).unwrap();
                if raw.len() < 16 || raw[0] != 2 { continue; }
                let mut nb = 0usize; let mut cb = 0usize; let mut bs = 0usize;
                unsafe { blosc_src::blosc_cbuffer_sizes(raw.as_ptr() as *const _, &mut nb, &mut cb, &mut bs) };
                let mut want = vec![0u8; nb];
                let got = unsafe {
                    blosc_src::blosc_decompress_ctx(raw.as_ptr() as *const _, want.as_mut_ptr() as *mut _, nb, 1)
                };
                assert_eq!(got as usize, nb, "{p}: {name}: c-blosc failed");
                assert_eq!(super::decompress(&raw).unwrap(), want, "{p}: {name}");
                n += 1;
            }
            eprintln!("{p}: {n} chunks identical to c-blosc");
        }
    }
}
