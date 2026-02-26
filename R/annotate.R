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
#' @param add_ref Additional reference(s) to add to the existing reference. Must
#'  be a list of SingleCellExperiment objects or something SingleR can handle.
#' @param add_labels Additional label(s) to add to the existing labels. Must be
#' a list of vectors.
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
                               add_ref,
                               add_labels,
                               ...) {
    if (species == "human") {
        ref_list <- make_human_celltype_ref_list(label_type)
    } else if (species == "mouse") {
        ref_list <- make_mouse_celltype_ref_list(label_type)
    }

    if ((missing(ref) || missing(labels)) && species == "") {
        stop("Please provide either ref/labels or species argument(s)")
    }

    if (!missing(add_ref) && !missing(add_labels)) {
        if (!is.list(add_ref) || !is.list(add_labels)) {
            stop("add_ref and add_labels must be lists")
        }

        ref_list$counts <- c(ref_list$counts, add_ref)
        ref_list$labels <- c(ref_list$labels, add_labels)
    }

    annotation <-
        SingleR::SingleR(test = Seurat::as.SingleCellExperiment(sobject),
                         ref = ref_list$counts,
                         labels = ref_list$labels,
                         aggr.ref = aggr_ref,
                         ...)

    sobject$cell_type <- annotation$labels

    # SingleR version 2.10.0 changed the cell_scores output from a matrix to a
    # DataFrame (not data.frame), so we need to check for both.
    # If it is a DataFrame, it will have a sub-DataFrame for each reference, so
    # we will take the max score across all references for each cell.
    if (is.matrix(annotation$scores)) {
        sobject$cell_scores <-
            apply(X = annotation$scores,
                  MARGIN = 1,
                  function(x) max(x, na.rm = TRUE))
    } else if (inherits(annotation$scores, "DFrame")) {
        sobject$cell_scores <-
            Reduce(cbind, annotation$scores) |>
                as.data.frame() |>
                dplyr::select(dplyr::starts_with("scores")) |>
                apply(1, max)
    }

    return(sobject)
}


##
#' Create human cell type reference list for annotation
#'
#' Builds a list of reference objects and their labels for human cell type
#' annotation using celldex references. Used internally by annotate_celltypes.
#'
#' @param label_type The label type to extract from each reference (e.g.,
#'   "label.main").
#' @return A list with two elements: 'counts' (list of reference objects) and
#'   'labels' (list of label vectors).
#' @noRd
make_human_celltype_ref_list <- function(label_type) {
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

    ref_info <-
        list(
            counts = list(
                hpca,
                huim,
                encode_bp
            ),
            labels = list(
                hpca[[label_type]],
                huim[[label_type]],
                encode_bp[[label_type]]
            )
        )

    return(ref_info)
}

##
#' Create mouse cell type reference list for annotation
#'
#' Builds a list of reference objects and their labels for mouse cell type
#' annotation using celldex references. Used internally by annotate_celltypes.
#'
#' @param label_type The label type to extract from each reference (e.g.,
#'   "label.main").
#' @return A list with two elements: 'counts' (list of reference objects) and
#'   'labels' (list of label vectors).
#' @noRd
make_mouse_celltype_ref_list <- function(label_type) {
    mord <- suppressMessages(celldex::MouseRNAseqData())
    moim <- suppressMessages(celldex::ImmGenData())
    moim$label.main <-
        stringr::str_remove_all(
            moim$label.main,
            c("B cells, pro" = "B cells", "DC" = "Dendritic cells")
        )

    ref_info <-
        list(
            counts = list(mord, moim),
            labels = list(mord[[label_type]], moim[[label_type]])
        )

    return(ref_info)
}