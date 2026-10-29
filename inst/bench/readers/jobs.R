# Profiling jobs: one reader call on one dataset each. A job is a function of
# the gDirSource backend; it runs the reader's own disk call.
# GI_DATA: root of the dataset folders (layout in README.md).
# GI_SCRATCH: where the CosMx fallback job keeps its decompressed matrix.
D <- path.expand(Sys.getenv("GI_DATA", "~/Documents/Datasets"))
scratch <- path.expand(Sys.getenv("GI_SCRATCH", tempdir()))
xen_lung <- file.path(D, "Xenium/LungCancer_10X/Xenium_V1_humanLung_Cancer_FFPE_outs")
xen_mini <- file.path(D, "Xenium/LungCancer_10X_small/Xenium_V1_humanLung_Cancer_FFPE_outs")
atera_bc <- file.path(D, "Atera/WTA_Preview_FFPE_Breast_Cancer_outs")
atera_kid <- file.path(D, "Atera/Atera_v1_WTA_Kidney_non_diseased_FFPE_outs")
cosmx_dir <- file.path(D, "CosMx/Breast_cancer/Flatfiles_RNA/flatFiles/BreastCancer")
vhd_mini <- file.path(D, "Cropped_datasets/VisiumHD/visiumhd_mini_crc/binned_outputs")
vhd_crc <- file.path(D, "Visium_HD/Human_Colorectal_Cancer")
st_mini <- file.path(D, "Cropped_datasets/Stereoseq/stereoseq_mini_mouse_eyeball")
st_testis <- file.path(D, "StereoSeq/C04687E314_backup/outs")

jobs <- list(
  xenium_mini_tx    = function(b) importXeniumDisk(xen_mini, backend = b)$load_transcripts(),
  xenium_mini_poly  = function(b) importXeniumDisk(xen_mini, backend = b)$load_polys(),
  xenium_mini_expr  = function(b) importXeniumDisk(xen_mini, backend = b)$load_expression(),
  xenium_lung_tx    = function(b) importXeniumDisk(xen_lung, backend = b)$load_transcripts(
                          path = file.path(xen_lung, "transcripts.parquet")),
  xenium_lung_tx_zarr = function(b) importXeniumDisk(xen_lung, backend = b)$load_transcripts(
                          path = file.path(xen_lung, "transcripts.zarr.zip")),
  xenium_lung_poly  = function(b) importXeniumDisk(xen_lung, backend = b)$load_polys(
                          path = file.path(xen_lung, "cell_boundaries.parquet")),
  xenium_lung_nuc   = function(b) importXeniumDisk(xen_lung, backend = b)$load_polys(
                          path = file.path(xen_lung, "nucleus_boundaries.parquet"), name = "nucleus",
                          part_col = NULL, split_geom = TRUE, split_geom_fmt = "nucleus_%d",
                          split_geom_sourcename = "cell_poly_id"),
  atera_bc_nuc      = function(b) importXeniumDisk(atera_bc, backend = b)$load_polys(
                          path = file.path(atera_bc, "nucleus_boundaries.parquet"), name = "nucleus",
                          part_col = NULL, split_geom = TRUE, split_geom_fmt = "nucleus_%d",
                          split_geom_sourcename = "cell_poly_id"),
  xenium_lung_expr  = function(b) importXeniumDisk(xen_lung, backend = b)$load_expression(),
  atera_bc_tx       = function(b) importXeniumDisk(atera_bc, backend = b)$load_transcripts(
                          path = file.path(atera_bc, "transcripts.parquet")),
  atera_bc_poly     = function(b) importXeniumDisk(atera_bc, backend = b)$load_polys(
                          path = file.path(atera_bc, "cell_boundaries.parquet")),
  atera_kid_expr    = function(b) importAteraDisk(atera_kid, backend = b)$load_expression(),
  cosmx_expr        = function(b) importCosMxDisk(cosmx_dir, backend = b)$load_expression(),
  # the R csv path: cosmxscan is only used for .gz input
  cosmx_expr_fallback = function(b) {
    plain <- file.path(scratch, "BreastCancer_exprMat_file.csv")
    if (!file.exists(plain)) {
      dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
      system2("gunzip", c("-c", shQuote(file.path(cosmx_dir, "BreastCancer_exprMat_file.csv.gz"))),
              stdout = plain)
    }
    importCosMxDisk(cosmx_dir, backend = b)$load_expression(path = plain)
  },
  cosmx_poly        = function(b) importCosMxDisk(cosmx_dir, backend = b)$load_polys(),
  vhd_mini_expr     = function(b) importVisiumHDDisk(vhd_mini, backend = b, bin = 2)$load_expression(),
  vhd_mini_tx       = function(b) importVisiumHDDisk(vhd_mini, backend = b, bin = 2)$load_transcripts(),
  vhd_expr          = function(b) importVisiumHDDisk(vhd_crc, backend = b, bin = 2)$load_expression(),
  vhd_tx            = function(b) importVisiumHDDisk(vhd_crc, backend = b, bin = 2)$load_transcripts(),
  st_mini_expr      = function(b) importStereoSeqDisk(st_mini, backend = b)$load_expression(),
  st_mini_tx        = function(b) importStereoSeqDisk(st_mini, backend = b)$load_binpoints(),
  st_expr           = function(b) importStereoSeqDisk(st_testis, backend = b)$load_expression(),
  st_tx             = function(b) importStereoSeqDisk(st_testis, backend = b)$load_binpoints()
)

# Inclusive-time stages, matched against Rprof call stacks
stages <- c(
  "tile plan (quadtree count scans)" = ".tile_plan_quadtree",
  "per-tile filtered scan (getTile)" = "getTile",
  "  collect into R"                 = "dplyr::collect",
  "R tile work (.store_write_points)" = ".store_write_points",
  "  data.table/data.frame copies"   = "copy",
  "  geometry frame (WKB etc.)"      = ".geom_parquet_frame",
  "  morton order"                   = ".morton_encode",
  "polygon build (terra/Giotto)"     = ".df_to_terra_poly",
  "parquet write (ExecPlan_Write)"   = "ExecPlan_Write",
  "zarr -> parquet"                  = ".zarr_ensure_parquet",
  "expression write loop"            = ".pestore_write_itr",
  "  write_parquet"                  = "write_parquet",
  "expression stats"                 = ".pestore_finalize_stats",
  "parquet intermediate (row_index)" = ".arrow_add_row_index",
  "fread"                            = "fread",
  "hdf5 reads"                       = "h5read",
  "forked children (wait)"           = "mclapply"
)

# MERSCOPE transcripts and Atera v1 zarr transcripts (the zarr input)
merscope_colon <- file.path(D, "Merscope/Colon_cancer/Vizgen_Colon_Cancer_V2")
atera_tiny <- file.path(D, "Atera/tiny_atera_dataset_outs")
jobs$merscope_tx <- function(b) importMerscopeDisk(merscope_colon, backend = b)$load_transcripts()
jobs$atera_tiny_tx <- function(b) importAteraDisk(atera_tiny, backend = b)$load_transcripts()
jobs$atera_kid_tx <- function(b) importAteraDisk(atera_kid, backend = b)$load_transcripts()
