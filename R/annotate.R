#' Annotate normal human and mouse cells in seurat object. Both human and mouse
#' cell reference are built by combining individual normal cell and immune
#' cell reference. This function is derived from the SingleR and hence can use
#' SingleR documentation for more help.
#'
#' @param sobject A seurat object that you would like to annotate
#' @param ref Pass "human" or "mouse" to annotate human or mouse cells
#'            respectively. You can pass your own ref if needed.
#' @param labels Automatically built in and passed as label.main column.
#'               You an also pass the metadata column that has the labels
#'               you are interested.
#' @param aggr.ref Arguments controlling the aggregation of the references
#'                 prior to annotation, see trainSingleR.
#' @param label_type Pick the label of interest from "label.main", "label.fine",
#'                   and "label.ont". The default is "label.main".
#' @param ... for other options in SingleR

annotate <- function(sobject,
                     ref,
                     labels,
                     aggr.ref,
                     label_type = "label.main",
                     ...) {
    if (ref == "human") {
        hpca <- celldex::HumanPrimaryCellAtlasData()
        huim <- celldex::MonacoImmuneData()
        hpca$label.main <-
            stringr::str_replace_all(hpca$label.main,
                                     c("T_cells" = "T cells",
                                       "B_cell" = "B cells",
                                       "NK_cell" = "NK cells",
                                       "Monocyte" = "Monocytes",
                                       "DC" = "Dendritic cells"))
        ref <- list(hpca,
                    huim)
        labels <- list(hpca[[label_type]],
                      huim[[label_type]])
    } else if (ref == "mouse") {
        mord <- celldex::MouseRNAseqData()
        moim <- celldex::ImmGenData()
        moim$label.main <-
            stringr::str_remove_all(moim$label.main,
                                    c("B cells, pro" = "B cells",
                                      "DC" = "Dendritic cells"))
        ref <- list(mord,
                    moim)
        labels <- list(mord[[label_type]],
                       moim[[label_type]])
    }
    annotation <-
        SingleR::SingleR(test = Seurat::as.SingleCellExperiment(object),
                         ref = ref,
                         labels = labels,
                         aggr.ref = aggr.ref,
                         ...)
    object$annotations <- annotation$labels
    object$cell_scores <-
        apply(X = annotation$scores,
              MARGIN = 1,
              function(x) max(x, na.rm = TRUE))
    return(object)
}