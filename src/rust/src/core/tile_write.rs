//! One point tile part as GiottoDisk's `.store_write_points()` writes it:
//! `row_index` (arrival order, 1-based), `x_index`, `y_index`, the WKB point
//! `geom`, then the metadata columns, rows in Morton order.

use std::collections::HashMap;
use std::fs::{create_dir_all, File};
use std::sync::Arc;

use arrow::array::{ArrayRef, BinaryBuilder, Float64Array, Int32Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::metadata::KeyValue;
use parquet::file::properties::WriterProperties;

/// The leading columns of a point tile part, rows taken in `ord`.
/// `row_offset` is added to the 1-based arrival position (parts after the
/// first continue the tile's numbering).
pub fn base_columns(
    x: &[f64], y: &[f64], ord: &[u32], row_offset: f64, with_geom: bool,
) -> (Vec<Field>, Vec<ArrayRef>) {
    let n = ord.len();
    let row_index = Int32Array::from_iter_values(ord.iter().map(|&i| (row_offset + i as f64 + 1.0) as i32));
    let xo = Float64Array::from_iter_values(ord.iter().map(|&i| x[i as usize]));
    let yo = Float64Array::from_iter_values(ord.iter().map(|&i| y[i as usize]));
    let mut fields = vec![
        Field::new("row_index", DataType::Int32, true),
        Field::new("x_index", DataType::Float64, true),
        Field::new("y_index", DataType::Float64, true),
    ];
    let mut cols: Vec<ArrayRef> = vec![Arc::new(row_index), Arc::new(xo), Arc::new(yo)];
    if with_geom {
        let mut geom = BinaryBuilder::with_capacity(n, n * 21);
        let mut wkb = [0u8; 21];
        wkb[0] = 1; // little endian
        wkb[1..5].copy_from_slice(&1u32.to_le_bytes()); // Point
        for &i in ord {
            wkb[5..13].copy_from_slice(&x[i as usize].to_le_bytes());
            wkb[13..21].copy_from_slice(&y[i as usize].to_le_bytes());
            geom.append_value(wkb);
        }
        fields.push(Field::new("geom", DataType::Binary, true));
        cols.push(Arc::new(geom.finish()));
    }
    (fields, cols)
}

/// Writes one parquet file (snappy, 1 Mi-row groups), creating its directory;
/// `geo` is the GeoParquet metadata, if any. It goes in the file's key-value
/// metadata, where GeoParquet readers (DuckDB spatial, SedonaDB) look for
/// it, as arrow C++ writes it; arrow-rs on its own keeps schema metadata
/// only inside the embedded `ARROW:schema`.
pub fn write_parquet(
    path: &str, fields: Vec<Field>, cols: Vec<ArrayRef>, geo: Option<String>,
) -> Result<(), String> {
    let md: HashMap<String, String> = geo.clone().map(|g| HashMap::from([("geo".to_string(), g)])).unwrap_or_default();
    let schema = Arc::new(Schema::new_with_metadata(fields, md));
    let batch = RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())?;
    if let Some(d) = std::path::Path::new(path).parent() {
        create_dir_all(d).map_err(|e| format!("{}: {e}", d.display()))?;
    }
    let f = File::create(path).map_err(|e| format!("{path}: {e}"))?;
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .set_max_row_group_size(1 << 20)
        .set_key_value_metadata(geo.map(|g| vec![KeyValue::new("geo".to_string(), g)]))
        .build();
    let mut w = ArrowWriter::try_new(f, schema, Some(props)).map_err(|e| e.to_string())?;
    w.write(&batch).map_err(|e| e.to_string())?;
    w.close().map_err(|e| e.to_string())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use parquet::file::reader::{FileReader, SerializedFileReader};

    /// GeoParquet readers find `geo` in the file's key-value metadata.
    #[test]
    fn geo_is_file_key_value_metadata() {
        let x = [1.0, 2.0];
        let y = [3.0, 4.0];
        let (fields, cols) = base_columns(&x, &y, &[0, 1], 0.0, true);
        let p = std::env::temp_dir().join(format!("gi_kv_{}.parquet", std::process::id()));
        let p = p.to_str().unwrap().to_string();
        write_parquet(&p, fields, cols, Some("{\"version\":\"1.0.0\"}".into())).unwrap();
        let r = SerializedFileReader::new(File::open(&p).unwrap()).unwrap();
        let kv = r.metadata().file_metadata().key_value_metadata().cloned().unwrap_or_default();
        std::fs::remove_file(&p).ok();
        let geo: Vec<_> = kv.iter().filter(|k| k.key == "geo").collect();
        assert_eq!(geo.len(), 1);
        assert_eq!(geo[0].value.as_deref(), Some("{\"version\":\"1.0.0\"}"));
    }
}
