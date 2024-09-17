test_that("tenx_load_qc pulls in a Seurat object", {
    sobj <- tenx_load_qc(test_path("test_10x_matrix"),
                         min_cells = 0,
                         min_features = 0,
                         violin_plot = FALSE)
    expect_true(inherits(sobj, "Seurat"))
})

test_that("process_seurat works", {
    sobj <- tenx_load_qc(test_path("test_10x_matrix"),
                         min_cells = 0,
                         min_features = 0,
                         violin_plot = FALSE)
    sobj <- process_seurat(sobj)
    expect_true(inherits(sobj, "Seurat"))
})

test_that("process_seurat works with resolution", {
    sobj <- tenx_load_qc(test_path("test_10x_matrix"),
                         min_cells = 0,
                         min_features = 0,
                         violin_plot = FALSE)
    sobj <- process_seurat(sobj, resolution = 0.8)
    expect_true(inherits(sobj, "Seurat"))
})
