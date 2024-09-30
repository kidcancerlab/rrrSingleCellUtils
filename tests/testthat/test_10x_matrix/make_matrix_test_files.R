# Making data to use for unit tests
library(Seurat)
library(rrrSingleCellUtils)
SeuratData::InstallData("hcabm40k")

test_p <-
    SeuratData::LoadData("hcabm40k") %>%
    UpdateSeuratObject() %>%
    process_seurat(resolution = 0.8) %>%
    PercentageFeatureSet(col.name = "percent_mito", pattern = "^MT-")

VlnPlot(test_p, features = c("percent_mito"), group.by = "orig.ident")

DimPlot(test_p)

# lets subset this down to one sample so we don't have duplicate barcodes once
# we cut the sample name off of the barcode
# Also, downsampling to reduce the size of the test data. I'm keeping the top
# variable genes so the data maintains some structure.
set.seed(1337)
cells_per_group <- 20
keep_n_genes <- 3000
test_p_sub <-
    subset(test_p, orig.ident == "MantonBM1") %>%
    subset(downsample = cells_per_group) %>%
    FindVariableFeatures(nfeatures = keep_n_genes)

mt_features <- grep("^MT-", rownames(test_p), value = TRUE)
test_p_sub <-
    test_p_sub[c(VariableFeatures(test_p_sub), mt_features), ] %>%
    process_seurat()

DimPlot(test_p_sub)

counts <- GetAssayData(test_p_sub, layer = "counts")
colnames(counts) <- stringr::str_remove(colnames(counts), ".+_[0-9]+-")

# Writing this out in two formats so we can test the read functions
DropletUtils::write10xCounts(counts,
                             path = "tests/testthat/test_10x_matrix/",
                             overwrite = TRUE)

DropletUtils::write10xCounts(counts,
                             path = "tests/testthat/test_10x_matrix/test.h5",
                             type = "HDF5")
