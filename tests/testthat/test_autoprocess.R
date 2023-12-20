test_that("process_sobj_gex works", {
    if(dir.exists("/home/gdrobertslab/lab/Counts/S0027/outs/filtered_feature_bc_matrix/")) {
        data <-
            tenx_load_qc("/home/gdrobertslab/lab/Counts/S0027/outs/filtered_feature_bc_matrix/",
                         violin_plot = FALSE,
                         mt_pattern = "^mm10-mt-|^hg19-MT-")

        sample_data <-
            readxl::read_xlsx("../sc_sequencing_log.xlsx") %>%
            dplyr::filter(Sample_ID == "S0027")
        sample_data$subset_nCount_RNA_min <- 10000
        sample_data$subset_nCount_RNA_max <- 75000
        sample_data$subset_percent.mt_max <- 20

        new_data <- process_sobj_gex(data, sample_data, "/gpfs0/scratch/mvc002")
        expect_true(class(new_data) == "Seurat")
    } else {
        skip("Required folder not available.")
    }
})
