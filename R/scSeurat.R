#' Load 10X data
#'
#' @param path10x Path to the 10X data to load
#' @param species Species data to analyze - see details
#' @param min_cells Passed to CreateSeuratObject : Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
#' @param min_features Passed to CreateSeuratObject : Include cells where at least this many features are detected.
#'
#' @details
#' For species:
#'    human - use for an alignment made to a human genome
#'    mouse - use for an alignment made to a murine genome
#'    mixHuman - use for an alignment made to a mixed genome, return human dataset
#'    mixMouse - use for an alignment made to a mixed genome, return murine dataset
#'
#' @return A \code{\link{Seurat}}
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- tenXLoadQC("path/to/10X/data/", species = "human")
#' }
tenXLoadQC <- function(path10x, species, min_cells = 5, min_features = 800) {
  raw.data <- Seurat::Read10X(path10x)

  if (species == "human") {
    mt_pattern <- "^MT-"
  } else if (species == "mouse") {
    mt_pattern <- "^mt-"
  } else if (species == "mixHuman") {
    mt_pattern <- "MT-"
    raw.data <- raw.data[grep(pattern = "^hg19", raw.data@Dimnames[[1]]), ]
    raw.data@Dimnames[[1]] <- substring(raw.data@Dimnames[[1]], 6)
  } else if (species == "mixMouse") {
    mt_pattern <- "^mt-"
    raw.data <- raw.data[grep(pattern = "^mm10", raw.data@Dimnames[[1]]), ]
    raw.data@Dimnames[[1]] <- substring(raw.data@Dimnames[[1]], 6)
  }

  seurat <- Seurat::CreateSeuratObject(raw.data,
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
