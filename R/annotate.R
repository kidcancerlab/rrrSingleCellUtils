#' Annotate normal human and mouse cells in seurat object.
#'
#' If specifying "species" argument, both human and mouse cell reference are
#' built by combining individual normal cell and immune celldex references.
#'
#' @param sobject A seurat object that you would like to annotate
#' @param species Takes "human" or "mouse" as a built in. Overrides ref and
#'  labels if set
#' @param ref Reference passed to Singler. Either a single or a list of multiple
#' @param labels Reference labels passed to SingleR. Either a single or a list
#'  of multiple
#' @param aggr_ref Passed to SingleR to control aggregation of the references
#'  prior to annotation, see trainSingleR.
#' @param label_type Either "label.main", "label.fine" or "label.ont" if using
#'  species argument. Not needed otherwise
#' @param ... Other options passed to SingleR
#'
#' @return A seurat object with cell_type and cell_scores added to metadata
#'
#' @export
annotate_celltypes <- function(sobject,
                               species = "",
                               ref,
                               labels,
                               aggr_ref = FALSE,
                               label_type = "label.main",
                               ...) {
    if (species == "human") {
        huim <- suppressMessages(celldex::MonacoImmuneData())
        huim$label.main <- stringr::str_replace_all(huim$label.main, "_", " ")

        hpca <- suppressMessages(celldex::HumanPrimaryCellAtlasData())
        hpca$label.main <-
            dplyr::case_match(
                hpca$label.main,
                .default = hpca$label.main,
                "Astrocyte"         ~ "Astrocytes",
                "B_cell"            ~ "B cells",
                "DC"                ~ "Dendritic cells",
                "NK_cell"           ~ "NK cells",
                "Macrophage"        ~ "Macrophages",
                "Monocyte"          ~ "Monocytes",
                "T_cells"           ~ "T cells"
            ) %>%
            stringr::str_replace_all("_", " ")

        encode_bp <- suppressMessages(celldex::BlueprintEncodeData())
        encode_bp$label.main <-
            dplyr::case_match(
                encode_bp$label.main,
                .default = encode_bp$label.main,
                "B-cells"           ~ "B cells",
                "CD4+ T-cells"      ~ "CD4+ T cells",
                "CD8+ T-cells"      ~ "CD8+ T cells",
                "DC"                ~ "Dendritic cells",
                "Endothelial cells" ~ "Endothelial_cells",
                "Epithelial cells"  ~ "Epithelial_cells",
                "Macrophagess"      ~ "Macrophages",
                "Smooth muscle"     ~ "Smooth_muscle_cells"
            ) %>%
            stringr::str_replace_all("_", " ")

        ref <- list(hpca,
                    huim,
                    encode_bp)
        labels <- list(hpca[[label_type]],
                       huim[[label_type]],
                       encode_bp[[label_type]])
    } else if (species == "mouse") {
        mord <- suppressMessages(celldex::MouseRNAseqData())
        moim <- suppressMessages(celldex::ImmGenData())
        moim$label.main <-
            stringr::str_remove_all(moim$label.main,
                                    c("B cells, pro" = "B cells",
                                      "DC" = "Dendritic cells"))
        ref <- list(mord,
                    moim)
        labels <- list(mord[[label_type]],
                       moim[[label_type]])
    }

    if (missing(ref) || missing(labels)) {
        stop("Please provide either ref/labels or species argument(s)")
    }

    annotation <-
        SingleR::SingleR(test = Seurat::as.SingleCellExperiment(sobject),
                         ref = ref,
                         labels = labels,
                         aggr.ref = aggr_ref,
                         ...)
    sobject$cell_type <- annotation$labels
    sobject$cell_scores <-
        apply(X = annotation$scores,
              MARGIN = 1,
              function(x) max(x, na.rm = TRUE))
    return(sobject)
}