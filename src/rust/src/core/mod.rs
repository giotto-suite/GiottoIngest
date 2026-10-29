//! What every input shares: GiottoDisk's per-tile row order, the point
//! tile part layout and writer, and the geoparquet metadata.

pub mod geo;
pub mod morton;
pub mod tile_write;
