#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param peak_beds DESCRIPTION.
#' @param min_peak_width DESCRIPTION.
#' @param max_peak_width DESCRIPTION.
#' @param frag_paths DESCRIPTION.
#' @param cell_ids DESCRIPTION.
#' @param n_regions_simul DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
merge_atac <- function(peak_beds,
                       min_peak_width = 20,
                       max_peak_width = 10000,
                       frag_paths,
                       cell_ids,
                       n_regions_simul = 2000) {
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
    obj_list <- lapply(seq_len(length(frag_paths)), make_chrom_seurat)
    seurat_obj <- merge(x = obj_list[[1]],
                        y = obj_list[-1],
                        add.cell.ids = cell_ids)

    return(seurat_obj)
}


#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param sobject DESCRIPTION.
#' @param gtf DESCRIPTION.
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


#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param sobject DESCRIPTION.
#' @param cutoff DESCRIPTION.
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



#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param sobject DESCRIPTION.
#' @param cutoff DESCRIPTION.
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


#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param sobject DESCRIPTION.
#' @param frag_file DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
calc_frip <- function(sobject, frag_file) {
    total_frag_df <- Signac::CountFragments(frag_file)
    
    sobject$total_frag <- total_frag_df$reads_count
    sobject$mononucleosomal <- total_frag_df$mononucleosomal
    sobject$nucleosome_free <- total_frag_df$nucleosome_free
    sobject <- Signac::FRiP(sobject,
                            assay = "ATAC",
                            total.fragments = "total_frag",
                            col.name = "FRiP")
    return(sobject)
}







# for sample in S0150 S0152 S0166 S0167 S0168 S0169 S0170
# do
#     echo ${sample}
#     python scripts/calcFragSize.py \
#         --fragFile ${base_path}/${sample}/atac_fragments.tsv.gz \
#         --quiet \
#         > output/atacFragLens/${sample}.txt
# done





# merged_datasets <- RNA_all[, colnames(RNA_all) %in% colnames(ATAC_all)] %>%
#     RunPCA()

# temp <- ATAC_all[, colnames(ATAC_all) %in% colnames(RNA_all)]
# merged_datasets[["ATAC"]] <- temp@assays$ATAC

# DefaultAssay(merged_datasets) <- "ATAC"

# merged_datasets <- RunTFIDF(merged_datasets) %>%
#     FindTopFeatures(min.cutoff = "q0") %>%
#     RunSVD() #%>%
#     #RunUMAP(reduction = "lsi", dims = 2:30)

# # Make joint UMAP
# merged_datasets2 <-
#     FindMultiModalNeighbors(
#     object = merged_datasets,
#     reduction.list = list("pca", "lsi"),
#     dims.list = list(1:50, 2:40),
#     verbose = TRUE)

# DefaultAssay(merged_datasets) <- "RNA"

# # build a joint UMAP visualization
# merged_datasets2 <- RunUMAP(object = merged_datasets2,
#                            nn.name = "weighted.nn",
#                            assay = "RNA") %>%
#     FindClusters(algorithm = 3, graph.name = "wsnn")