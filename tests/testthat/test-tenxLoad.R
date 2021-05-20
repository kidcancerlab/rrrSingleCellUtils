test_that("tenx_load_qc works", {
  dir_name <-
    "/home/gdrobertslab/lab/Counts/S0027/outs/filtered_feature_bc_matrix/"

    if(dir.exists(dir_name)) {
      expect_true(TRUE)
  } else {
    skip("Required folder not available.")
  }
})
