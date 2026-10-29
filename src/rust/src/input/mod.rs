//! Inputs, one module per on-disk format, feeding the shared core.
//!
//! - `parquet`: Xenium / Atera preview / MERSCOPE transcripts, bin points
//!   and polygon vertices; rows routed to a GiottoDisk quadtree plan.
//! - `zarr`: Atera v1 `transcripts.zarr.zip`, decoded per grid cell into
//!   GiottoDisk's Atera tile plan.
//! - `csv` (next): delimited and wide matrix inputs.

pub mod parquet;
pub mod zarr;
