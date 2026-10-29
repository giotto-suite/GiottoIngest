test_that("the build reports its capabilities", {
    expect_true(has_capability("points_parquet"))
    expect_true(has_capability("points_atera_zarr"))
    expect_false(has_capability("expression_csv"))
    expect_false(has_capability(NA_character_))
    expect_type(gi_capabilities(), "character")
})
