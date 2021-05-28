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
#'     species (eg, "^mm10" or "^hg19")
#' @param violin_plot If TRUE, produces a violin plot
#'
#' @return A \code{\link{Seurat}}
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- tenXLoadQC("path/to/10X/data/", species_pattern = "^mm9")
#' }
tenx_load_qc <- function(path_10x, min_cells = 5, min_features = 800,
                         mt_pattern = "^mt-|^MT-", species_pattern = "",
                         violin_plot = TRUE) {
  raw_data <- Seurat::Read10X(path_10x)

  if (species_pattern != "") {
    raw_data <- raw_data[grep(pattern = species_pattern,
                              raw_data@Dimnames[[1]]), ]
    raw_data@Dimnames[[1]] <- substring(raw_data@Dimnames[[1]], 6)
  }

  seurat <- Seurat::CreateSeuratObject(raw_data,
                                       min.cells = min_cells,
                                       min.features = min_features)
  seurat <- Seurat::PercentageFeatureSet(seurat,
                                         pattern = mt_pattern,
                                         col.name = "percent.mt")

  if (violin_plot == TRUE) {
    print(Seurat::VlnPlot(seurat,
                          features = c("nFeature_RNA",
                                       "nCount_RNA",
                                       "percent.mt"),
                          ncol = 3))
  }

  return(seurat)
}


#' Extract cellecta barcode information from a sam or bam file
#'
#' @param file Input sam or bam file
#' @param verbose Optional updates during parsing
#' @param output File to write out
#' @param samtools_module For slurm, which modules to load
#'
#' @return A tibble of the cell ids and lineage tracing barcodes
#' @export
#'
#' @details For some reason, this function does not cancel when ordered. Either
#' kill the R session or kill the process directly on the server
#'
#' @examples
#' \dontrun{
#' cid_lt <- gen_cellecta_bc_data(file = "path/to/file.bam",
#'                                verbose = TRUE,
#'                                samtools_module = "GCC/9.3.0 SAMtools/1.10")
#'
#' cid_lt <- gen_cellecta_bc_data(
#'    file = "/home/gdrobertslab/lab/Counts/S0027/outs/S0016-S0027-bc2.sam",
#'    verbose = TRUE,
#'    samtools_module = "GCC/9.3.0 SAMtools/1.10")
#' }
gen_cellecta_bc_data <- function(file, verbose = FALSE, output = tempfile(),
                                 samtools_module = FALSE) {

  package_dir <- find.package("rrrSingleCellUtils")

  system_cmd <- paste(package_dir,
                      "/exec/getCellectaBarcodes.pl --sam ",
                      file,
                      " > ", output,
                      sep = "")

  if (verbose) {
    system_cmd <- paste(system_cmd, "-v")
  }

  if (samtools_module != FALSE) {
    system_cmd <- paste("ml load ", samtools_module, "; ",
                        system_cmd, "; ",
                        "ml unload ", samtools_module,
                        sep = "")
  } else if (Sys.which("samtools") == "") {
    stop("\"samtools\" command not found in $PATH, do you need to provide module
            info for slurm?")
  }

  system(system_cmd)

  results <- readr::read_delim(output, delim = "\t", col_names = TRUE)

  return(results)
}

#' Process Lineage Tracing Barcodes
#'
#' @param sobject Name of the Seurat object that contains the cell barcodes
#'     (cids) that will be matched and integrated
#' @param histogram Will trigger the function to generate and output a
#'     histogram plot of the top 40 most frequent lineage tracing barcodes
#' @param ymax Upper limit of the y axis (ie, for creating side-by-side
#'     comparisons)
#' @param relative This will normalize cell counts to total number of cells
#'     containing barcodes
#' @param cid_lt Table of cellecta cell id and lineage tracing barcodes.
#'     Returned by gen_cellecta_bc_data()
#' @param col_fill Fill color
#' @param ret_list If TRUE, return a list of barcode frequencies instead of a
#'     Seurat object
#' @param title Title of the histogram, if generated
#'
#' @return either a \code{\link{Seurat}} object or a data frame of barcode
#'     frequences
#' @export
#'
#' @examples
#' \dontrun{
#' cid_lt <- gen_cellecta_bc_data(file = "path/to/file.bam",
#'                                verbose = TRUE,
#'                                samtools_module = "GCC/9.3.0 SAMtools/1.10")
#' output <- process_ltbc(sobject, cid_lt = cid_lt, histogram = TRUE)
#' }
process_ltbc <- function(sobject, cid_lt, histogram = FALSE,
                         col_fill = "#4DBBD5FF", ymax = NA, relative = FALSE,
                         title = "LT Barcode Frequency", ret_list = FALSE) {

  # Deduplicate redundant reads
  cid_lt <- cid_lt %>%
    dplyr::distinct() %>%

    # Match extracted barcode reads against the Cellecta barcode tables
    dplyr::mutate(label14 = bc14f[lt_14],
                  label30 = bc30f[lt_30]) %>%
    dplyr::mutate(label14 = stringr::str_remove(label14, ".+-"),
                  label30 = stringr::str_remove(label30, ".+-")) %>%

    # Eliminate barcodes that don't match the Cellecta barcode tables
    stats::na.omit() %>%

    # Concatenate the two barcodes into a single compound column
    dplyr::mutate(label = paste(label14, label30, sep = "-")) %>%
    dplyr::pull(label, cid) %>%
    sort()

  # Integrate the lineage tracing barcode into the Seurat object metadata
  sobject$lt <- cid_lt[sobject@assays$RNA@counts@Dimnames[[2]]]

  # Generate the frequency tables
  ylabel <- "Number of Cells"
  bc_freq <- as.data.frame(table(sobject$lt)) %>%
    dplyr::rename(freq = Freq) %>%
    dplyr::arrange(-freq)

  if (isTRUE(relative)) {
    bc_freq$freq <- bc_freq$freq / sum(bc_freq$freq) * 100
    ylabel <- "Percentage of Cells"
  }

  bc_plot_data <- utils::head(bc_freq, n = 40)

  # Create histogram graphs (default using blue color from npg from ggsci)
  if (isTRUE(histogram)) {
    print(ggplot2::ggplot(bc_plot_data, ggplot2::aes(x = stats::reorder(Var1,
                                                                        -freq),
                                                     y = freq)) +
            ggplot2::geom_bar(fill = col.fill, stat = "identity") +
            ggplot2::ggtitle(title) +
            ggplot2::ylab(ylabel) +
            ggplot2::xlab("Lineage Barcode") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                               hjust = 1)) +
            ggplot2::ylim(0, ymax))
  }

  if (isTRUE(ret_list)) {
    return(bc_freq)
  } else {
    return(sobject)
  }
}

#' Regress out cell cycle effects
#'
#' @param sobject Seurat object to be processed
#' @param cc_regress If set to Y, the process with run without user input and
#'     will automatically proceed to cell cycle regression
#' @param find_pcs Number of principal components to generate in the re-do PCA
#'     post-CC regression
#' @param use_pcs Number of principal components to use in the post-regression
#'     dimensional reduction
#' @param use_res Resolution to input to FindClusters
#' @param method Type of dimensional reduction to use, currently supports either
#'     umap or tsne
#'
#' @return A Seurat object
#' @export
#'
#' @details
#' The kill_cc function will identify cell cycle components within a dataset.
#' After an initial scoring using the Seurat CellCycleScoring function, the user
#' will be shown a dimensional reduction plot with cells labeled by cell cycle.
#' If indicated, the user can then trigger a process to regress out the effects
#' of cell cycle within the dataset.  The function will then proceed to re-do
#' the PCA and jackstraw if needed, then show a dimensional reduction plot
#' post-regression and retun the corrected Seurat object.
#' Input must be a Seurat object that already has PCA and dimensional reduction
#' data (umap or tsne) attached.
#'
#'
#' @examples
#' \dontrun{
#' load("~/analyses/roberts/dev/rrrSingleCellUtils/testData/test_cc.RData")
#' test <- kill_cc(os, use_pcs = 5, cc_regress = "Y")
#' }
kill_cc <- function(sobject, cc_regress = "N", find_pcs = 20, use_pcs = 3,
                    use_res = 0.5, method = "umap") {
  sobject <- Seurat::CellCycleScoring(sobject,
                                      s.features = Seurat::cc.genes$s.genes,
                                      g2m.features = Seurat::cc.genes$g2m.genes,
                                      set.ident = TRUE)

  if (method != "umap" & method != "tsne") {
    print("You need to set a method supported by this function")
  }

  plot_cc <- Seurat::DimPlot(sobject,
                             reduction = method,
                             label = TRUE,
                             pt.size = 1)

  print(plot_cc)

  if (cc_regress != "Y") {
    cc_regress <-
      readline(prompt =
                 "Proceed with regression of cell cycle-dependent genes (Y/N)?")
  }

  if (cc_regress == "Y") {
    sobject <- Seurat::ScaleData(sobject,
                                 vars.to.regress = c("S.Score", "G2M.Score"),
                                 features = rownames(x = sobject))
    sobject <- Seurat::RunPCA(sobject, npcs = find_pcs)

    if (method == "tsne") {
      sobject <- Seurat::RunTSNE(sobject, reduction = "pca", dims = 1:use_pcs)
    } else if (method == "umap") {
      sobject <- Seurat::RunUMAP(sobject, reduction = "pca", dims = 1:use_pcs)
    }

    sobject <- Seurat::FindNeighbors(sobject,
                                     reduction = "pca",
                                     dims = 1:use_pcs)

    sobject <- Seurat::FindClusters(sobject, resolution = use_res)

    print(Seurat::DimPlot(sobject,
                          reduction = method,
                          label = TRUE,
                          pt.size = 1,
                          group.by = "Phase"))

    print(Seurat::DimPlot(sobject,
                          reduction = method,
                          label = TRUE,
                          pt.size = 1))

  } else {
    print("No CC regression performed.")
  }
  return(sobject)
}

#' Plot Cell Cycle Distribution
#'
#' @param sobject A Seurat object
#' @param plot_type Type of plot to output, one of "bar" or "pie"
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_cc(marrow)
#'
#' plot_cc(marrow) + facet_wrap(~ Cluster, nrow = 1)
#' }
plot_cc <- function(sobject, plot_type = "bar") {
  if (!"Phase" %in% names(sobject@meta.data)) {
    stop("Phase info not available in this object.
         Run Seurat::CellCycleScoring() to estimate cell cycle phases")
  }

  cid2 <- tibble::tibble(Cluster = sobject$seurat_clusters,
                         CellId = names(sobject$seurat_clusters)) %>%
    dplyr::full_join(tibble::tibble(Phase = sobject$Phase,
                                    CellId = names(sobject$Phase)),
                     by = "CellId") %>%
    dplyr::group_by(Phase, Cluster) %>%
    dplyr::summarize(Freq = length(Phase), .groups = "drop") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Cluster) %>%
    dplyr::mutate(Proportion = Freq / sum(Freq)) %>%
    dplyr::ungroup()

  plot_obj <- ggplot2::ggplot(cid2,
                              ggplot2::aes(x = Cluster,
                                           y = Proportion,
                                           fill = Phase)) +
    ggplot2::geom_bar(width = 1, stat = "identity", color = "white") +
    ggplot2::scale_fill_brewer(palette = "Blues") +
    ggplot2::theme_minimal()

  if (plot_type == "pie") {
    plot_obj <- ggplot2::ggplot(cid2,
                                ggplot2::aes(x = "",
                                             y = Proportion,
                                             fill = Phase)) +
      ggplot2::geom_bar(stat = "identity", color = "white") +
      ggplot2::coord_polar("y", start = 0) +
      ggplot2::scale_fill_brewer(palette = "Blues") +
      ggplot2::theme_minimal() +
      Seurat::NoAxes() +
      ggplot2::facet_wrap(~ Cluster)
  }
  return(plot_obj)
}
