#' Merge multiple ATAC-seq samples into a single Seurat object
#'
#' FUNCTION_DESCRIPTION
#'
#' @param peak_beds DESCRIPTION.
#' @param min_peak_width DESCRIPTION.
#' @param max_peak_width DESCRIPTION.
#' @param frag_paths DESCRIPTION.
#' @param cell_ids DESCRIPTION.
#' @param n_regions_simul DESCRIPTION.
#' @param threads Number of threads to use.
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
merge_atac <- function(peak_beds,
                       min_peak_width = 20,
                       max_peak_width = 10000,
                       frag_paths,
                       cell_ids,
                       n_regions_simul = 2000,
                       threads = 1) {
    message("Make sure that paths and cell_ids are in the same order")
    message("peak_beds: ", paste(peak_beds, sep = " "))
    message("frag_paths: ", paste(frag_paths, sep = " "))
    message("cell_ids: ", paste(cell_ids, sep = " "))

    # read peaks into a dataframe and reduce them
    reduced_peaks <- readr::read_tsv(peak_beds,
                                     comment = "#",
                                     col_names = c("chr", "start", "end"),
                                     show_col_types = FALSE) %>%
        GenomicRanges::makeGRangesFromDataFrame() %>%
        Signac::reduce()

    # filter out peaks that are too small or too large
    reduced_peaks <- reduced_peaks[
        BiocGenerics::width(reduced_peaks) >= min_peak_width &
        BiocGenerics::width(reduced_peaks) <= max_peak_width
                                ]

    # Function to get fragment files, count and make Seurat object
    make_chrom_seurat <- function(i) {
        frags <- Signac::CreateFragmentObject(path = frag_paths[[i]])

        counts <- Signac::FeatureMatrix(fragments = frags,
                                        features = reduced_peaks,
                                        process_n = n_regions_simul)

        obj <- Signac::CreateChromatinAssay(counts = counts,
                                            fragments = frags) %>%
                    SeuratObject::CreateSeuratObject(assay = "ATAC")
    }

    # Use apply to make the Seurat objects and then merge them
    obj_list <-
        parallel::mclapply(seq_len(length(frag_paths)),
                           make_chrom_seurat,
                           mc.cores = threads)

    seurat_obj <- merge(x = obj_list[[1]],
                        y = obj_list[-1],
                        add.cell.ids = cell_ids)

    return(seurat_obj)
}


#' Add annotation to ATAC data in a Seurat object
#'
#' FUNCTION_DESCRIPTION
#'
#' @param sobject Seurat object to be processed
#' @param gtf String of path to a gtf file.
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
annotate_atac <- function(sobject, gtf) {
    if (is.null(gtf)) {
        message("No gtf file provided")
        stop
    }

    # Create annotation to put into seurat object
    annotations <-
        GenomicFeatures::makeTxDbFromGFF(gtf,
                                         format = "gtf") %>%
        GenomicFeatures::transcripts(columns = c("tx_id", "tx_name"))

    package_dir <- find.package("rrrSingleCellUtils")

    gene_info <-
        system(paste0("python ",
                      package_dir,
                      "/exec/parseGtfForAnnot.py --gtf ",
                      gtf),
               intern = TRUE) %>%
        read.table(text = ., header = TRUE) %>%
        dplyr::full_join(tibble::tibble(tx_name = annotations$tx_name,
                                        tx_id = annotations$tx_id)) %>%
        dplyr::arrange(tx_id)

    # "This should all be TRUE"
    if (all(gene_info$tx_name == annotations$tx_name)) {

        annotations$gene_biotype <- gene_info$gene_biotype
        annotations$gene_id <- gene_info$gene_id
        annotations$gene_name <- gene_info$gene_name

        Signac::Annotation(sobject) <- annotations

        return(sobject)
    } else {
        message("parsing gtf file failed")
        stop()
    }
}


#' Add nucleosome signal to ATAC data in a Seurat object
#'
#' FUNCTION_DESCRIPTION
#'
#' @param sobject Seurat object to be processed
#' @param cutoff DESCRIPTION.
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
add_nucleosome_signal <- function(sobject, cutoff = 4) {
    sobject <- Signac::NucleosomeSignal(sobject)
    sobject$nucleosome_group <- ifelse(sobject$nucleosome_signal > cutoff,
                                       paste0("NS > ", cutoff),
                                       paste0("NS < ", cutoff))
    return(sobject)
}

#' Calculate TSS enrichment for a Seurat object
#'
#' FUNCTION_DESCRIPTION
#'
#' @param sobject Seurat object to be processed
#' @param cutoff DESCRIPTION.
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
tss_enrichment <- function(sobject, cutoff = 2) {
    sobject <- Signac::TSSEnrichment(sobject, fast = FALSE)
    sobject$high_tss <- ifelse(sobject$TSS.enrichment > cutoff,
                               "High",
                               "Low")
    return(sobject)
}

#' Calculate FRiP for a Seurat object
#'
#' FUNCTION_DESCRIPTION
#'
#' @param sobject Seurat object to be processed
#' @param frag_files Named list of fragment files. The names should be the
#'    string to be prepended to the cell barcodes.
#' @param verbose Should functions be verbose?
#'
#' @export
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
calc_frip <- function(sobject, frag_files, verbose = FALSE) {
    # Get fragments for each sample in order and add sample name to CB column
    # if only one file provided, with no name
    if (is.null(names(frag_files)) & length(frag_files) == 1) {
        message("Did you mean to add a name to this list?")
        total_frag_df <- Signac::CountFragments(frag_files[[1]],
                                                verbose = verbose)
    # if only one file provided, with a name, append name to CB column
    } else if (!is.null(names(frag_files)) & length(frag_files) == 1) {
        total_frag_df <-
            Signac::CountFragments(frag_files[[1]],
                                   verbose = verbose) %>%
            dplyr::mutate(CB = paste(names(frag_files), CB, sep = "_"))
    # Otherwise, if names are provided, append to CB column
    } else if (!is.null(names(frag_files))) {
        total_frag_df <- data.frame(CB = character())
        for (frag in names(frag_files)) {
            if (verbose) message("\nProcessing ", frag, " fragments")
            total_frag_df <-
                Signac::CountFragments(frag_files[[frag]],
                                       verbose = verbose) %>%
                dplyr::mutate(CB = paste(frag, CB, sep = "_")) %>%
                dplyr::bind_rows(total_frag_df)
        }
    # If no names and more than one file, throw error
    } else {
        stop("Please provide a named list of fragment files")
    }

    # Ensure that total_frag_df is in the same order as the cells in sobject
    total_frag_df <-
        total_frag_df %>%
        dplyr::filter(CB %in% colnames(sobject)) %>%
        dplyr::arrange(match(CB, colnames(sobject)))

    if (nrow(total_frag_df) != ncol(sobject)) {
        warning("Number of cells in sobject and total_frag_df do not match")
        warning("Check that your names are correct in your frag_files list. ",
                "Note that this function adds the underscore before the cell ",
                "barcode.")
        return(total_frag_df)
    }

    # Populate sobject with metadata
    sobject$total_frag <- total_frag_df$reads_count
    sobject$mononucleosomal <- total_frag_df$mononucleosomal
    sobject$nucleosome_free <- total_frag_df$nucleosome_free
    sobject <- Signac::FRiP(sobject,
                            assay = "ATAC",
                            total.fragments = "total_frag",
                            col.name = "FRiP",
                            verbose = verbose)
    return(sobject)
}


#' Process a Seurat object with ATAC data
#'
#' FUNCTION_DESCRIPTION
#'
#' @param sobject Seurat object to be processed
#' @param assay Assay that has the ATAC data
#' @param verbose Should functions be verbose?
#' @param umap_dims PCA dimensions to use for UMAP
#'
#' @export
#'
#' @return A processed Seurat object
#' @examples
#' # ADD_EXAMPLES_HERE
process_seurat_atac <- function(sobject,
                                assay = "ATAC",
                                verbose = FALSE,
                                umap_dims = 2:30) {
    # Save current assay so we can reset it later
    old_active_ident <- Seurat::DefaultAssay(sobject)

    sobject <-
        sobject %>%
        Signac::RunTFIDF(verbose = verbose) %>%
        Signac::FindTopFeatures(min.cutoff = "q0", verbose = verbose) %>%
        Signac::RunSVD(verbose = verbose) %>%
        Seurat::FindNeighbors(reduction = "lsi", verbose = verbose) %>%
        Seurat::FindClusters(algorithm = 3, verbose = verbose) %>%
        Seurat::RunUMAP(reduction = "lsi",
                        dims = umap_dims,
                        verbose = verbose)

    # Reset active assay
    Seurat::DefaultAssay(sobject) <- old_active_ident
    return(sobject)
}

# atac_frag_size <- function(frag_files, verbose = FALSE) {
#     sys_cmd <- paste0("python scripts/calcFragSize.py --fragFile ",
#                       )
#     if (verbose) {
        
#     }
# }
# for sample in S0150 S0152 S0166 S0167 S0168 S0169 S0170
# do
#     echo ${sample}
#     python scripts/calcFragSize.py \
#         --fragFile ${base_path}/${sample}/atac_fragments.tsv.gz \
#         --quiet \
#         > output/atacFragLens/${sample}.txt
# done

#' Merge GEX and ATAC data
#'
#' FUNCTION_DESCRIPTION
#'
#' @param gex_sobj Seurat object with GEX data to be processed
#' @param atac_sobj Seurat object with ATAC data to be processed
#' @param atac_assay_name Name of the assay in atac_sobj that has the ATAC data
#' @param gex_pca_dims GEX PCA dimensions to use for UMAP
#' @param atac_pca_dims ATAC PCA dimensions to use for UMAP
#' @param verbose Should functions be verbose?
#'
#' @export
#'
#' @return A processed Seurat object
#' @examples
#' # ADD_EXAMPLES_HERE
merge_gex_atac <- function(gex_sobj,
                           atac_sobj,
                           atac_assay_name = "ATAC",
                           gex_pca_dims = 1:dim(Seurat::Embeddings(gex_sobj))[2],
                           atac_pca_dims = 1:dim(Seurat::Embeddings(atac_sobj))[2],
                           verbose = FALSE) {
    # Make seurat object from GEX keeping only cells present in ATAC
    merged_data <-
        gex_sobj[, colnames(gex_sobj) %in% colnames(atac_sobj)] %>%
        Seurat::RunPCA()

    # Make temporary seurat object from ATAC with cells present in GEX
    temp <-
        atac_sobj[, colnames(atac_sobj) %in% colnames(gex_sobj)] %>%
        Signac::RunTFIDF(verbose = verbose) %>%
        Signac::FindTopFeatures(min.cutoff = "q0",
                                verbose = verbose) %>%
        Signac::RunSVD(verbose = verbose)

    # Stuff it into the merged_data object
    merged_data[["ATAC"]] <- temp@assays[[atac_assay_name]]

    # Make joint UMAP
    merged_data <-
        merged_data %>%
        Seurat::FindMultiModalNeighbors(reduction.list = list("pca", "lsi"),
                                        dims.list = list(gex_pca_dims, atac_pca_dims),
                                        verbose = verbose)

    # build a joint UMAP visualization
    merged_data <-
        merged_data %>%
        Seurat::RunUMAP(nn.name = "weighted.nn",
                        assay = "RNA",
                        verbose = verbose) %>%
        Seurat::FindClusters(algorithm = 3,
                             graph.name = "wsnn",
                             verbose = verbose)

    return(merged_data)
}