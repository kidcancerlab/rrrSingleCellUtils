test_that("tenx_load_qc works", {
  if(dir.exists(
    "/home/gdrobertslab/lab/Counts/S0027/outs/filtered_feature_bc_matrix/")) {
    expect_true(is)
  } else {
    skip("Required folder not available.")
  }
})
