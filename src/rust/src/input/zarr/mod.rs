//! Zarr v2 inputs: the zip archive reader and the Atera v1 transcript
//! decoder (a port of GiottoDisk's `.atera_tx_cell_worker()`).

pub mod archive;
pub mod atera;
mod atera_layout;
