library(Seurat)

# Test case 1: FDL on default graph with weighted edges
pbmc_small <- run_fdl(pbmc_small)
expect_true("fdl" %in% Seurat::Reductions(pbmc_small))

# Test case 2: FDL on non-existent graph
tryCatch({
  pbmc_small <- run_fdl(pbmc_small, graph = "non_existent_graph")
}, error = function(e) {
  expect_true(grepl("non_existent_graph graph not found", e$message))
})
