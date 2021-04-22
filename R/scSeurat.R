#' Load 10X data
#'
#' @param path_10x Path to 10X RNAseq data "filtered_feature_bc_matrix" folder
#' @param min_cells Passed to CreateSeuratObject: Include features detected in
#'     at least this many cells. Will subset the counts matrix as well. To
#'     reintroduce excluded features, create a new object with a lower cutoff.
#' @param min_features Passed to CreateSeuratObject: Include cells where at
#'     least this many features are detected.
#' @param mt_pattern Pattern used to identify mitochondrial reads
#' @param species_pattern Pattern used to select only reads from a single
#'     species
#'
#' @return A \code{\link{Seurat}}
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- tenXLoadQC("path/to/10X/data/", species_pattern = "^mm9")
#' }
tenx_load_qc <- function(path_10x, min_cells = 5, min_features = 800,
                         mt_pattern = "^mt-|^MT-", species_pattern = "^hg19") {
  raw_data <- Seurat::Read10X(path_10x)
  raw_data <- raw_data[grep(pattern = species_pattern,
                            raw_data@Dimnames[[1]]), ]
  raw_data@Dimnames[[1]] <- substring(raw_data@Dimnames[[1]], 6)

  seurat <- Seurat::CreateSeuratObject(raw_data,
                               min.cells = min_cells,
                               min.features = min_features)
  seurat <- Seurat::PercentageFeatureSet(seurat,
                                 pattern = mt_pattern,
                                 col.name = "percent.mt")

  print(Seurat::VlnPlot(seurat,
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3))

  return(seurat)
}
