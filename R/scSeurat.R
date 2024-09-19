#' Load 10X data
#'
#' @param path_10x Path to 10X RNAseq data "filtered_feature_bc_matrix" folder
#' @param h5_file Path to 10X h5 file
#' @param frag_file Path to 10X ATAC fragments file. Only required for ATAC data
#' @param min_cells Passed to CreateSeuratObject: Include features detected in
#'     at least this many cells. Will subset the counts matrix as well. To
#'     reintroduce excluded features, create a new object with a lower cutoff.
#' @param min_features Passed to CreateSeuratObject: Include cells where at
#'     least this many features are detected.
#' @param mt_pattern Pattern used to identify mitochondrial reads (eg, "^MT-"
#'     Must add species_pattern if remove_species_pattern = FALSE
#'     (eg, "^hg19-MT-" or "^hg19-MT-|^mm10-mt-")
#' @param species_pattern Pattern used to select only reads from a single
#'     species (eg, "^mm10-" or "^hg19-")
#' @param exp_type Experiment type (one of "GEX", "ATAC" or "GEX+ATAC")
#' @param violin_plot If TRUE (default), produces a violin plot
#' @param remove_species_pattern Specifies if want to remove species_pattern
#'     prefix from gene names. If TRUE (default), removes species_pattern
#'     prefix.
#' @param sample_name Sample name. Used in violin plots for title
#'
#' @return A \code{\link{Seurat}}
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- tenXLoadQC("path/to/10X/data/", species_pattern = "^mm9-")
#' }
tenx_load_qc <- function(path_10x = "",
                         h5_file = "",
                         frag_file = stringr::str_replace(h5_file,
                                                          "filtered_feature_bc_matrix.h5",
                                                          "atac_fragments.tsv.gz"),
                         min_cells = 5,
                         min_features = 800,
                         mt_pattern = "^mt-|^MT-",
                         species_pattern = "",
                         exp_type = "GEX",
                         remove_species_pattern = TRUE,
                         violin_plot = TRUE,
                         sample_name = path_10x) {

    check_load_inputs(remove_species_pattern,
                      species_pattern,
                      mt_pattern,
                      path_10x,
                      h5_file,
                      exp_type)

    # Load in the data using either the 5h file or the 10x folder
    if (h5_file != "") {
        raw_data <- suppressWarnings(Seurat::Read10X_h5(h5_file))
        if ("Peaks" %in% names(raw_data)) {
            atac_raw_data <- raw_data$Peaks
        }
        if ("Gene Expression" %in% names(raw_data)) {
            rna_raw_data <- raw_data$`Gene Expression`
        } else {
            rna_raw_data <- raw_data
        }
    } else {
        rna_raw_data <- Seurat::Read10X(path_10x)
        if ("Gene Expression" %in% names(rna_raw_data)) {
            rna_raw_data <- rna_raw_data[["Gene Expression"]]
        }
    }

    if (exp_type == "GEX+ATAC" && min_features != 0) {
        min_features <- 0
        message("Setting min_features to 0 to load both GEX and ATAC data.\n",
                "Carefully filter your data downstream.")
    }

    if (grepl("GEX", exp_type)) {
        gex_orig_cells <- nrow(rna_raw_data)
        gex_first_ten_genes <- rownames(rna_raw_data)[1:10]
        # subset the data to only include the species of interest
        rna_raw_data <-
            filter_raw_data(rna_raw_data,
                            species_pattern,
                            remove_species_pattern)

        if (nrow(rna_raw_data) == 0) {
            stop("No genes left in object after filtering using ",
                 "species_pattern! Check species_pattern argument. ",
                 "Before filtering gex data using species_pattern the first ",
                 "genes were: ",
                 paste(gex_first_ten_genes, collapse = "\n"))
        }

        seurat <-
            Seurat::CreateSeuratObject(rna_raw_data,
                                       min.cells = min_cells,
                                       min.features = min_features)

        if (ncol(seurat) == 0) {
            stop("No data left after applying min.cells and min.features.",
                 " Check your data and arguments.")
        }

        seurat <-
            Seurat::PercentageFeatureSet(seurat,
                                         pattern = gsub("_", "-", mt_pattern),
                                         col.name = "percent.mt")
            # Need to change all underscores to dashes due to CreateSeuratObject
            # doing the same

        if (sum(seurat$percent.mt, na.rm = TRUE) == 0) {
            warning("No mitochondrial reads found!")
            warning("If you have a sample aligned to a mixed reference, make ",
                    "sure that your species_pattern and mt_pattern arguments ",
                    "are appropriate.")
            print(paste("Potential mitochondrial genes:",
                        rownames(seurat)[grep(rownames(seurat),
                                         pattern = "mt",
                                         ignore.case = TRUE)]))
        }
    }

    if (grepl("ATAC", exp_type)) {
        # subset the data to only include the species of interest
        atac_first_ten_peaks <- head(rownames(rna_raw_data), 10)

        atac_raw_data <-
            filter_raw_data(atac_raw_data,
                            species_pattern,
                            remove_species_pattern)

        if (nrow(atac_raw_data) == 0) {
            stop("No peaks left in object after filtering using ",
                 "species_pattern! Check species_pattern argument. ",
                 "Before filtering gex data using species_pattern the first ",
                 "peaks were: ",
                 paste(atac_first_ten_peaks, collapse = "\n"))
        }


        if (exp_type == "ATAC") {
            seurat <-
                Signac::CreateChromatinAssay(counts = atac_raw_data,
                                             sep = c(":", "-"),
                                             fragments = frag_file,
                                             min.cells = min_cells) %>%
                Seurat::CreateSeuratObject(assay = "ATAC")
        } else if (exp_type == "GEX+ATAC") {
            seurat[["ATAC"]] <-
                Signac::CreateChromatinAssay(counts = atac_raw_data,
                                             sep = c(":", "-"),
                                             fragments = frag_file,
                                             min.cells = min_cells)
        } else {
            stop("This shouldn't be possible. I have no idea how you got here.")
        }

        seurat <-
            Seurat::PercentageFeatureSet(seurat,
                                         pattern = gsub("_", "-", mt_pattern),
                                         col.name = "percent_mt_atac",
                                         assay = "ATAC")

    }



    if (violin_plot && grepl("GEX", exp_type)) {
        print(Seurat::VlnPlot(seurat,
                            features = c("nFeature_RNA",
                                         "nCount_RNA",
                                         "percent.mt"),
                            ncol = 3) +
            patchwork::plot_annotation(title = path_10x,
                                       theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 25))))
    }

    return(seurat)
}

check_load_inputs <- function(remove_species_pattern,
                              species_pattern,
                              mt_pattern,
                              path_10x,
                              h5_file,
                              exp_type) {

    # If removing species names and species pattern is in mt_pattern
    if (remove_species_pattern == TRUE &&
        grepl(species_pattern %>% stringr::str_remove_all("\\^"),
              mt_pattern) &&
        species_pattern != "") {
        warning("\nDon't put species prefix in mt_pattern when
                remove_species_pattern == TRUE\n",
                immediate. = TRUE)
        stop()
    }

    # If not removing species names and species pattern not in mt_pattern
    if (remove_species_pattern == FALSE &&
        grepl(species_pattern %>% stringr::str_remove_all("\\^"),
            mt_pattern) == FALSE) {
        warning("\nMake sure your mt_pattern have species prefixes in it when
                remove_species_pattern == FALSE
                mt_pattern should look like: ^hg19-MT-\n",
                immediate. = TRUE)
        stop()
    }

    if (!exp_type %in% c("GEX", "ATAC", "GEX+ATAC")) {
        warning("\nexp_type must be one of 'GEX', 'ATAC', or 'GEX+ATAC'\n",
                immediate. = TRUE)
        stop()
    }

    if (grepl("ATAC", exp_type) && h5_file == "") {
        warning("\nYou need to specify a h5_file for ATAC data\n",
                immediate. = TRUE)
        stop()
    }

    if (path_10x != "" && h5_file != "") {
        warning("\nYou can't specify both a path_10x and h5_file\n",
                immediate. = TRUE)
        stop()
    }

    if (h5_file != "") {
        if (suppressWarnings(try(system("which h5pfc",
                                        intern = TRUE,
                                        ignore.stderr = TRUE))) %>%
                length() == 0) {
            warning("\nYou need to have h5pfc installed to use the h5_file ",
                    "argument. Perhaps ml load HDF5 before you start R?\n",
                    immediate. = TRUE)
            stop()
        }
    }
}

filter_raw_data <- function(raw_data_matrix,
                            species_pattern,
                            remove_species_pattern) {
    if (species_pattern != "") {
        raw_data_matrix <-
            raw_data_matrix[grep(pattern = gsub("-", "_", species_pattern),
                            rownames(raw_data_matrix)), ]
        # Then remove the species_pattern prefix from the peak names
        if (remove_species_pattern) {
            rownames(raw_data_matrix) <-
                sub(gsub("-",
                            "_",
                            species_pattern),
                    "",
                    rownames(raw_data_matrix))
        }
    }
    return(raw_data_matrix)
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

  results <-
    readr::read_delim(output,
                      delim = "\t",
                      col_names = TRUE,
                      show_col_types = FALSE)

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
  sobject$lt <-
    cid_lt[stringr::str_remove(Seurat::Cells(sobject),
                               "-1$")] %>%
    as.vector()


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
            ggplot2::geom_bar(fill = col_fill, stat = "identity") +
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
#'     will automatically proceed to cell cycle regression. If set to Ask, will
#'     prompt the user. If set to N no regression will be performed.
#' @param show_plots Should the plots be printed?
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
kill_cc <- function(sobject,
                    cc_regress = "N",
                    show_plots = TRUE,
                    find_pcs = 20,
                    use_pcs = 3,
                    use_res = 0.5,
                    method = "umap") {
  sobject <- Seurat::CellCycleScoring(sobject,
                                      s.features = Seurat::cc.genes$s.genes,
                                      g2m.features = Seurat::cc.genes$g2m.genes,
                                      set.ident = TRUE)

  if (method != "umap" & method != "tsne") {
    print("You need to set a method supported by this function")
  }

  plot_cc <-
    Seurat::DimPlot(sobject,
                    reduction = method,
                    label = TRUE,
                    pt.size = 1) +
    theme_roberts()

  if (show_plots) {
    print(plot_cc)
  }

  if (cc_regress == "Ask") {
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

    if (show_plots) {
        print(Seurat::DimPlot(sobject,
                            reduction = method,
                            label = TRUE,
                            pt.size = 1,
                            group.by = "Phase") +
                theme_roberts())

        print(Seurat::DimPlot(sobject,
                            reduction = method,
                            label = TRUE,
                            pt.size = 1) +
                theme_roberts())
    }

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

#' Process a Seurat object
#'
#' @inheritParams Seurat::RunUMAP
#' @inheritParams Seurat::FindClusters
#' @param sobject A Seurat object
#' @param run_umap_dims Number of PCA dimensions to use in RunUMAP()
#'      and FindNeighbors()
#' @param graph_name Name of graph to use for the clustering algorithm in FindClusters()
#' @param neighbor_k_param Number of neighbors (k.param) to use in FindNeighbors()
#' @param umap_n_neighbors Number of neighbors (n.neighbors) to use in RunUMAP()
#' @param umap_metric Metric (metric) to use in RunUMAP()
#'
#' @return A Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' test <- process_seurat(pbmc_small)
#' }
process_seurat <- function(sobject,
                           verbose = FALSE,
                           run_umap_dims = 1:30,
                           assay = "RNA",
                           resolution = 0.3,
                           reduction = "pca",
                           graph_name = "RNA_snn",
                           neighbor_k_param = 30,
                           umap_n_neighbors = 30L,
                           umap_metric = "euclidean") {

    sobject <-
        sobject %>%
        Seurat::NormalizeData(verbose = verbose,
                              assay = assay) %>%
        Seurat::FindVariableFeatures(verbose = verbose,
                                     assay = assay) %>%
        Seurat::ScaleData(verbose = verbose,
                          assay = assay) %>%
        Seurat::RunPCA(verbose = verbose,
                       assay = assay)

    # Make sure run_umap_dims isn't more than we have in pca
    if (max(run_umap_dims) > ncol(Seurat::Embeddings(sobject,
                                                     reduction = reduction))) {
        run_umap_dims <- seq_len(ncol(Seurat::Embeddings(sobject,
                                                 reduction = reduction)))
        warning(paste("run_umap_dims argument is greater than the number of",
                      "dimensions in the reduction. Using all dimensions."))
    }

    sobject <-
        sobject %>%
        Seurat::RunUMAP(dims = run_umap_dims,
                        reduction = reduction,
                        n.neighbors = umap_n_neighbors,
                        metric = umap_metric,
                        verbose = verbose,
                        assay = assay) %>%
        Seurat::FindNeighbors(dims = run_umap_dims,
                              reduction = reduction,
                              verbose = verbose,
                              assay = assay,
                              k.param = neighbor_k_param) %>%
        Seurat::FindClusters(resolution = resolution,
                             verbose = verbose,
                             graph.name = graph_name)
}

#' Use median and standard deviation to subset a Seurat object based on specific features
#'
#' @param sobject Seurat object
#' @param sd_down Number of standard deviations below the median to subset
#' @param sd_up Number of standard deviations above the median to subset
#' @param make_plots Whether to make plots of the features before and after subsetting
#' @param features Vector of features to use for subsetting
#' @param sample_name Name of sample to use in plot titles
#'
#' @return A Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' filtered <- auto_subset(SeuratObject::pbmc_small)
#' }
auto_subset <- function(sobject,
                        sd_down = 1,
                        sd_up = 2,
                        make_plots = TRUE,
                        features = c("nCount_RNA",
                                     "nFeature_RNA"),
                        sample_name = NULL) {

    cutoff_table <-
        sobject@meta.data %>%
        dplyr::select(dplyr::all_of(features)) %>%
        tidyr::pivot_longer(cols = dplyr::everything(),
                            names_to = "feature",
                            values_to = "value") %>%
        dplyr::group_by(feature) %>%
        dplyr::summarize(median_val = stats::median(value),
                         sd_val = stats::sd(value),
                         .groups = "drop") %>%
        dplyr::mutate(min_val = median_val - sd_down * sd_val,
                      max_val = median_val + sd_up * sd_val)

    if (make_plots) {
        print(feature_hist(sobject,
                           features = features,
                           cutoff_table = cutoff_table))
    }

    for (feature_name in features) {
        cutoffs <-
            cutoff_table %>%
            dplyr::filter(feature == feature_name)

        values <- Seurat::FetchData(object = sobject, vars = feature_name)
        sobject <-
            sobject[, which(x = values >= cutoffs$min_val[1] &
                                values <= cutoffs$max_val[1])]
    }

    # if (make_plots) {
    #     print(feature_hist(sobject,
    #                        features = features,
    #                        cutoff_table = cutoff_table))
    # }
    return(sobject)
}
