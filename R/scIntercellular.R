#' Create a gene list containing putative targets of ligand activity.
#'
#' @details This gene list is needed as one of the inputs for the
#'    ligand-receptor analysis. This list should be generated to identify genes
#'    in the same target cells used in the findLigands function as the "targets"
#'    parameter. You do not have to use this particular function to generate
#'    that gene list, but you can use this function to do it.
#' @param sobject A Seurat object containing all of the cells for analysis
#'    (required)
#' @param id1 The idents of the target cells (ligand-stimulated cells)
#'    (required)
#' @param id2 The idents of the unstimlated cells for comparison (required)
#' @param pval The p-value to use as a cutoff for up-regulation of genes
#'    (default = 0.05)
#' @param logfc The log2 fold-change to use as cutoff for up-regulation of
#'    genes (default = 0.25)
#' @param spec The species of the gene set (default = "human", can also be
#'    "mouse")
#' @export
#'
#' @return A list of target genes
#'
#' @examples
#' \dontrun{
#' targets <- find_tar_genes(seurat_obj,
#'                           id1 = "d14",
#'                           id2 = "d35",
#'                           logfc = 0.25,
#'                           spec = "human")
#' }

find_tar_genes <- function(sobject, id1, id2, pval = 0.05, logfc = 0.25,
                         spec = "human") {

  if (spec == "human") {
    ligand_list <- human_ligand_list
  } else if (spec == "mouse") {
    ligand_list <- mouse_ligand_list
  } else {
    stop("for findTarGenes, spec must be defined as either human or mouse")
  }

  targets <-
    Seurat::FindMarkers(sobject, ident.1 = id1, ident.2 = id2) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(p_val_adj <= pval & abs(avg_log2FC) >= logfc) %>%
    dplyr::pull(gene) %>%
    dplyr::intersect(ligand_list)

  return(targets)
}
