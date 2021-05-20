
## find_ligands()
test_that("find_ligands() works", {
  if(file.exists("testData/find_ligands_test.RData")) {
    load("testData/find_ligands_test.RData")

    output <- find_ligands(combo,
                           gset = targets_2,
                           receiver = "Anchors",
                           senders = c("M1-like",
                                       "M2-like",
                                       "AlvMac",
                                       "DC"),
                           gset_spec = "human",
                           rec_spec = "human",
                           send_spec = "mouse")

  } else {
    skip("Required data not available.")
  }
})


load("testData/find_tar_genes_test.RData")

## find_tar_genes()
targets <- find_tar_genes(os,
                          id1 = "d14",
                          id2 = "d35",
                          logfc = 0.25,
                          spec = "human")

test_that("find_tar_genes output is a character vector", {
  expect_true(is.character(targets))
})

test_that("find_tar_genes output is a the proper length", {
  expect_length(object = targets, n = 2399)
})
