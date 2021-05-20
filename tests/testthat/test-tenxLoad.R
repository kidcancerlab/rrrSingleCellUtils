test_that("tenx_load_qc works", {
  dir_name <-
    "/home/gdrobertslab/lab/Counts/S0027/outs/filtered_feature_bc_matrix/"

  if(dir.exists(dir_name)) {
    seurat_obj <- tenx_load_qc(dir_name, violin_plot = FALSE)
    
    test_that("tenx_load_qc loads in data", {
      expect_true(class(seurat_obj) == "Seurat")
    })
    
  } else {
    skip("Required folder not available.")
  }
})
