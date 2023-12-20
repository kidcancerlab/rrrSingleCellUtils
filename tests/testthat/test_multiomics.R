s0150_dir <- "/home/gdrobertslab/lab/Counts/S0150/"
h_arc_gtf <- "/home/gdrobertslab/lab/GenRef/10x-human_arc/genes/genes.gtf"


test_that("annotate_atac works", {
    if (dir.exists(paste0(s0150_dir, "filtered_feature_bc_matrix/"))) {
        data <-
            suppressWarnings(tenx_load_qc(h5_file = paste0(s0150_dir,
                                                          "filtered_feature_bc_matrix.h5"),
                                          frag_file = paste0(s0150_dir,
                                                             "atac_fragments.tsv.gz"),
                                          species_pattern = "^hg19-",
                                          exp_type = "GEX+ATAC",
                                          violin_plot = FALSE))

        Seurat::DefaultAssay(data) <- "ATAC"

        new_data <-
            annotate_atac(data,
                          gtf = h_arc_gtf)

        expect_true(class(new_data) == "Seurat")
    } else {
        skip("Required folder not available.")
    }
})


test_that("add_atac_metadata works", {
    if (dir.exists(paste0(s0150_dir, "filtered_feature_bc_matrix/"))) {
        data <-
            suppressWarnings(tenx_load_qc(h5_file = paste0(s0150_dir,
                                                          "filtered_feature_bc_matrix.h5"),
                                          frag_file = paste0(s0150_dir,
                                                             "atac_fragments.tsv.gz"),
                                          species_pattern = "^hg19-",
                                          exp_type = "GEX+ATAC",
                                          violin_plot = FALSE))

        new_data <-
            add_atac_metadata(data,
                              gtf = h_arc_gtf,
                              frag_files = paste0(s0150_dir,
                                                  "atac_fragments.tsv.gz"))

        expect_true(class(new_data) == "Seurat")
    } else {
        skip("Required folder not available.")
    }
})

# Test calc_frip
test_that("calc_frip works", {
    if (dir.exists(paste0(s0150_dir, "filtered_feature_bc_matrix/"))) {
        data <-
            suppressWarnings(tenx_load_qc(h5_file = paste0(s0150_dir,
                                                          "filtered_feature_bc_matrix.h5"),
                                          frag_file = paste0(s0150_dir,
                                                             "atac_fragments.tsv.gz"),
                                          species_pattern = "^hg19-",
                                          exp_type = "GEX+ATAC",
                                          violin_plot = FALSE))

        Seurat::DefaultAssay(data) <- "ATAC"

        new_data <-
            annotate_atac(data,
                          gtf = h_arc_gtf)

        new_data <-
            calc_frip(new_data)

        expect_true(class(new_data) == "Seurat")
    } else {
        skip("Required folder not available.")
    }
})
