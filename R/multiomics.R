#' Merge multiple ATAC-seq samples into a single Seurat object
#'
#' Reads peak BED files, reduces overlapping peaks, and creates a merged
#' Seurat object with ATAC-seq data from multiple samples.
#'
#' @param peak_beds Character vector of paths to BED files containing peak
#'     regions for each sample.
#' @param min_peak_width Minimum width (in base pairs) for peaks to be
#'     included in the analysis. Default is 20.
#' @param max_peak_width Maximum width (in base pairs) for peaks to be
#'     included in the analysis. Default is 10000.
#' @param frag_paths Character vector of paths to fragment files (one per
#'     sample).
#' @param cell_ids Character vector of sample identifiers to prepend to cell
#'     barcodes. Must be in the same order as frag_paths.
#' @param n_regions_simul Number of regions to process simultaneously when
#'     building the feature matrix. Default is 2000.
#' @param threads Number of threads to use for parallel processing. Default
#'     is 1.
#'
#' @export
#'
#' @return A Seurat object with a merged ATAC assay containing data from all
#'     samples.
#' @examples
#' \dontrun{
#' merged_obj <- merge_atac(
#'   peak_beds = c("sample1_peaks.bed", "sample2_peaks.bed"),
#'   frag_paths = c("sample1_fragments.tsv.gz", "sample2_fragments.tsv.gz"),
#'   cell_ids = c("sample1", "sample2")
#' )
#' }
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
        reduced_peaks@ranges@width >= min_peak_width &
        reduced_peaks@ranges@width <= max_peak_width
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

        return(obj)
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
#' Creates genomic annotations from a GTF file and adds them to the ATAC
#' assay in a Seurat object. The function extracts transcript information,
#' gene biotypes, and gene identifiers from the GTF file.
#'
#' @param sobject Seurat object containing ATAC data to be annotated. The
#'     default assay should be set to "ATAC" or "peaks" before calling this
#'     function.
#' @param gtf Character string specifying the path to a GTF file containing
#'     genomic annotations.
#'
#' @export
#'
#' @return A Seurat object with genomic annotations added to the ATAC assay.
#' @examples
#' \dontrun{
#' atac_obj <- annotate_atac(
#'   sobject = atac_seurat,
#'   gtf = "path/to/genes.gtf"
#' )
#' }
annotate_atac <- function(sobject,
                          gtf) {
    if (is.null(gtf)) {
        message("No gtf file provided")
        stop
    }

    if (!Seurat::DefaultAssay(sobject) %in% c("ATAC", "peaks")) {
        warning("If your default assay isn't your ATAC data, this function",
                "will fail. Make sure to set DefaultAssay() properly.")
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
        dplyr::arrange(tx_id) %>%
        dplyr::filter(tx_name %in% annotations@elementMetadata$tx_name)
        # the filter accounts for any transcripts that are dropped while
        # creating the "annotations" object

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
#' Calculates nucleosome signal for each cell using Signac and adds it as
#' metadata. Also creates a categorical grouping variable based on the
#' specified cutoff.
#'
#' @param sobject Seurat object containing ATAC data to be processed.
#' @param cutoff Numeric value for the nucleosome signal cutoff. Cells with
#'     nucleosome signal above this value are grouped as "NS > cutoff" and
#'     cells below are grouped as "NS < cutoff". Default is 4.
#'
#' @export
#'
#' @return A Seurat object with nucleosome signal added as metadata in the
#'     nucleosome_signal and nucleosome_group columns.
#' @examples
#' \dontrun{
#' atac_obj <- add_nucleosome_signal(sobject = atac_seurat, cutoff = 4)
#' }
add_nucleosome_signal <- function(sobject,
                                  cutoff = 4) {
    sobject <- Signac::NucleosomeSignal(sobject)
    sobject$nucleosome_group <- ifelse(sobject$nucleosome_signal > cutoff,
                                       paste0("NS > ", cutoff),
                                       paste0("NS < ", cutoff))
    return(sobject)
}

#' Calculate TSS enrichment for a Seurat object
#'
#' Calculates transcription start site (TSS) enrichment scores for each cell
#' using Signac. Also creates a categorical variable indicating high or low
#' TSS enrichment based on the specified cutoff.
#'
#' @param sobject Seurat object containing ATAC data to be processed.
#' @param cutoff Numeric value for the TSS enrichment cutoff. Cells with TSS
#'     enrichment above this value are classified as "High" and cells below
#'     are classified as "Low". Default is 2.
#'
#' @export
#'
#' @return A Seurat object with TSS enrichment added as metadata in the
#'     TSS.enrichment and high_tss columns.
#' @examples
#' \dontrun{
#' atac_obj <- tss_enrichment(sobject = atac_seurat, cutoff = 2)
#' }
tss_enrichment <- function(sobject,
                           cutoff = 2) {
    sobject <- Signac::TSSEnrichment(sobject, fast = FALSE)
    sobject$high_tss <- ifelse(sobject$TSS.enrichment > cutoff,
                               "High",
                               "Low")
    return(sobject)
}

#' Calculate FRiP for a Seurat object
#'
#' Calculates the fraction of reads in peaks (FRiP) metric for each cell in
#' an ATAC-seq Seurat object. This function counts fragments from the
#' provided fragment files and calculates quality metrics including total
#' fragments, mononucleosomal fragments, nucleosome-free fragments, and FRiP.
#'
#' @param sobject Seurat object containing ATAC data to be processed.
#' @param frag_files Named list of paths to fragment files. The names should
#'     correspond to sample identifiers that will be prepended to cell
#'     barcodes (with an underscore separator). If providing a single
#'     unnamed fragment file, cell names must follow standard CellRanger
#'     format (e.g., "ATGC-1").
#' @param verbose Logical indicating whether functions should produce
#'     detailed output messages. Default is FALSE.
#'
#' @export
#'
#' @return A Seurat object with FRiP and fragment count metadata added. If
#'     there is a mismatch between cells in the object and fragment files, a
#'     data frame of fragment counts is returned instead for inspection.
#' @examples
#' \dontrun{
#' atac_obj <- calc_frip(
#'   sobject = atac_seurat,
#'   frag_files = list(sample1 = "path/to/fragments1.tsv.gz",
#'                     sample2 = "path/to/fragments2.tsv.gz")
#' )
#' }
calc_frip <- function(sobject,
                      frag_files,
                      verbose = FALSE) {
    # Get fragments for each sample in order and add sample name to CB column
    # if only one file provided, with no name for the list and the cell names
    # don't look like normal cellranger names, throw a warning
    if (all(stringr::str_detect(colnames(sobject),
                                        "^[ATGC]+-[0-9]+$"))) {
        # all good, I think
        total_frag_df <- Signac::CountFragments(frag_files[[1]],
                                                verbose = verbose)
    } else if (is.null(names(frag_files)) &&
                length(frag_files) == 1) {
        message("Did you mean to add a name to this list? Your cell names ",
                "don't look like normal cellranger names.")
        stop()
    # if only one file provided, with a name, append name to CB column
    } else if (!is.null(names(frag_files)) && length(frag_files) == 1) {
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

#' Add ATAC specific metadata to the Seurat object
#'
#' Add nucleosome signal, TSS enrichment, and FRiP to a Seurat object
#'
#' @param sobject Seurat object to be processed
#' @param gtf String of path to a gtf file.
#' @param nucl_cutoff Cutoff for nucleosome signal
#' @param tss_cutoff Cutoff for TSS enrichment
#' @param frag_files Named list of paths to fragment files.
#' @param verbose Should functions be verbose?
#'
#' @export
#'
#' @return A Seurat object
add_atac_metadata <- function(sobject,
                              gtf,
                              nucl_cutoff = 4,
                              tss_cutoff = 2,
                              frag_files,
                              verbose = TRUE) {
    sobject <-
        annotate_atac(sobject,
                      gtf = gtf) %>%
        add_nucleosome_signal(cutoff = nucl_cutoff) %>%
        tss_enrichment(cutoff = tss_cutoff) %>%
        calc_frip(frag_files = frag_files,
                  verbose = verbose)
}

#' Process a Seurat object with ATAC data
#'
#' Performs standard ATAC-seq processing workflow including TF-IDF
#' normalization, feature selection, dimensionality reduction using singular
#' value decomposition (SVD/LSI), clustering, and UMAP generation.
#'
#' @param sobject Seurat object containing ATAC data to be processed.
#' @param assay Name of the assay containing the ATAC data. Default is
#'     "ATAC".
#' @param verbose Logical indicating whether processing functions should
#'     produce detailed output messages. Default is FALSE.
#' @param umap_dims Integer vector specifying which LSI dimensions to use
#'     for UMAP calculation. Default is 2:30 (dimensions 2 through 30).
#' @param resolution Numeric value for clustering resolution passed to
#'     FindClusters. Higher values result in more clusters. Default is 0.3.
#' @param reduction Character string specifying which dimensionality
#'     reduction to use for clustering and UMAP. Default is "lsi".
#'
#' @export
#'
#' @return A processed Seurat object with LSI reduction, UMAP coordinates
#'     (stored as "umap_atac"), and cluster assignments. The original active
#'     assay is restored after processing.
#' @examples
#' \dontrun{
#' processed_atac <- process_seurat_atac(
#'   sobject = atac_seurat,
#'   resolution = 0.5
#' )
#' }
process_seurat_atac <- function(sobject,
                                assay = "ATAC",
                                verbose = FALSE,
                                umap_dims = 2:30,
                                resolution = 0.3,
                                reduction = "lsi") {
    # Save current assay so we can reset it later
    old_active_ident <- Seurat::DefaultAssay(sobject)

    sobject <-
        sobject %>%
        Signac::RunTFIDF(verbose = verbose) %>%
        Signac::FindTopFeatures(min.cutoff = "q0", verbose = verbose) %>%
        Signac::RunSVD(verbose = verbose) %>%
        Seurat::FindNeighbors(reduction = reduction, verbose = verbose) %>%
        Seurat::FindClusters(algorithm = 3,
                             verbose = verbose,
                             resolution = resolution) %>%
        Seurat::RunUMAP(reduction = reduction,
                        dims = umap_dims,
                        verbose = verbose,
                        reduction.name = "umap_atac")

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
#' Combines gene expression (GEX) and ATAC-seq data from separate Seurat
#' objects into a single multimodal object. Performs joint dimensionality
#' reduction and clustering using weighted nearest neighbor (WNN) analysis.
#'
#' @param gex_sobj Seurat object containing gene expression data. Should
#'     have PCA reduction already computed.
#' @param atac_sobj Seurat object containing ATAC-seq data. Should have LSI
#'     reduction already computed.
#' @param atac_assay_name Character string specifying the name of the ATAC
#'     assay in atac_sobj. Default is "ATAC".
#' @param gex_pca_dims Integer vector specifying which PCA dimensions to use
#'     from the GEX data. Default uses all available PCA dimensions.
#' @param atac_pca_dims Integer vector specifying which LSI dimensions to
#'     use from the ATAC data. Default uses all available LSI dimensions.
#' @param verbose Logical indicating whether processing functions should
#'     produce detailed output messages. Default is FALSE.
#'
#' @export
#'
#' @return A Seurat object containing both RNA and ATAC assays with joint
#'     UMAP coordinates and WNN-based cluster assignments. Only cells present
#'     in both input objects are retained.
#' @examples
#' \dontrun{
#' multimodal_obj <- merge_gex_atac(
#'   gex_sobj = gex_seurat,
#'   atac_sobj = atac_seurat,
#'   gex_pca_dims = 1:30,
#'   atac_pca_dims = 2:30
#' )
#' }
merge_gex_atac <- function(gex_sobj,
                           atac_sobj,
                           atac_assay_name = "ATAC",
                           gex_pca_dims =
                               1:dim(Seurat::Embeddings(gex_sobj))[2],
                           atac_pca_dims =
                               1:dim(Seurat::Embeddings(atac_sobj,
                                                        reduction = "lsi"))[2],
                           verbose = FALSE) {
    # Make seurat object from GEX keeping only cells present in ATAC
    merged_data <-
        gex_sobj[, colnames(gex_sobj) %in% colnames(atac_sobj)] %>%
        Seurat::RunPCA()

    # Make temporary seurat object from ATAC with cells present in GEX
    temp <-
        atac_sobj[, colnames(atac_sobj) %in% colnames(gex_sobj)]

    # Stuff it into the merged_data object and process it
    merged_data[["ATAC"]] <- temp@assays[[atac_assay_name]]
    Seurat::DefaultAssay(merged_data) <- "ATAC"
    merged_data <-
        merged_data %>%
        Signac::RunTFIDF(verbose = verbose) %>%
        Signac::FindTopFeatures(min.cutoff = "q0",
                                verbose = verbose) %>%
        Signac::RunSVD(verbose = verbose)

    # Make joint UMAP
    merged_data <-
        merged_data %>%
        Seurat::FindMultiModalNeighbors(
            reduction.list = list("pca", "lsi"),
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
