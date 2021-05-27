usethis::use_data(bc14f,
                  bc30f,
                  human_ligand_list,
                  mouse_ligand_list,
                  internal = TRUE,
                  compress = "xz")

# From http://satijalab.org/seurat/cell_cycle_vignette.html
exp.mat <- read.table(file = "nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, as.is = TRUE, row.names =1)

marrow <- CreateSeuratObject(counts = exp.mat) %>%
  NormalizeData(.) %>%
  FindVariableFeatures(., selection.method = "vst") %>%
  ScaleData(., features = rownames(.)) %>%
  RunPCA(., features = VariableFeatures(.)) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.3) %>%
  RunUMAP(dims = 1:20)

save(marrow, file = "testData/kill_cc_test.RData", compress = "xz")
