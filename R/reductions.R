#' Run Force-Directed Layout (FDL) on a Seurat object
#'
#' This function runs Force-Directed Layout (FDL) on a Seurat object.
#' FDL is a graph layout algorithm that positions nodes based on attractive and
#'  repulsive forces between them.
#'
#' @param sobject The Seurat object to run FDL on
#' @param graph The name of the graph in the Seurat object to use for FDL
#' @param weighted Logical indicating whether the graph edges should be weighted
#'
#' @export
#'
#' @return The modified Seurat object with the FDL results added
#'
#' @examples
#' \dontrun{
#' # Load Seurat object
#' data("pbmc_small")
#'
#' # Run FDL on the Seurat object
#' pbmc_small <- run_fdl(pbmc_small)
#'
#' # Access the FDL results
#' fdl_results <- pbmc_small[["fdl"]]
#' }
run_fdl <- function(sobject,
                    graph = "RNA_snn",
                    weighted = TRUE) {
    if (!graph %in% names(sobject@graphs)) {
        stop(graph, " graph not found in Seurat object")
    }

    graph <-
        igraph::graph_from_adjacency_matrix(
            adjmatrix = sobject[[graph]],
            mode = "undirected",
            weighted = weighted,
            add.colnames = TRUE
        )

    fdl <- igraph::layout_with_fr(graph, grid = "nogrid")
    rownames(fdl) <- colnames(sobject)
    colnames(fdl) <- c("fdl_1", "fdl_2")

    sobject[["fdl"]] <-
        SeuratObject::CreateDimReducObject(
            embeddings = fdl,
            key = "fdl_",
            assay = SeuratObject::DefaultAssay(sobject)
        )

    return(sobject)
}