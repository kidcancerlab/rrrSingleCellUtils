# two_sobj <- qs::qread("/home/gdrobertslab/mvc002/analyses/roberts/dev/testSnps/output/patient_met_1/two_sobj.qs")

# c_b_t <-
#     two_sobj@meta.data %>%
#     select(used_clusters) %>%
#     rename(cell_group = used_clusters) %>%
#     rownames_to_column("cell_barcode") %>%
#     as_tibble() %>%
#     mutate(bam_file = paste0("/home/gdrobertslab/lab/Counts/",
#                         str_remove(cell_barcode, "_.*"),
#                         "/outs/possorted_genome_bam.bam"),
#             cell_barcode = str_remove(cell_barcode, ".+_"))


#' Generate a tree using SNPs derived from single cell clusters
#'
#' @param cellid_bam_table A tibble with three columns: cell_barcode,
#'  cell_group and bam_file. The cell_barcode column should contain the cell
#'  barcode, the cell_group column should contain the cluster label and the
#'  bam_file column should contain the path to the bam file for that cell.
#' @param temp_dir The directory to write temporary files to.
#' @param slurm_base The directory to write slurm output files to.
#' @param account The hpc account to use.
#' @param ploidy Either a file path to a ploidy file, or a string indicating
#'  the ploidy.
#' @param ref_fasta The path to the reference fasta file.
#' @param min_depth The minimum depth to use when calling SNPs.
#' @param submit Whether to submit the sbatch jobs to the cluster or not.
#'
#' @return A hclust tree
#' @export
#'
#' @examples
#' \dontrun{
#' placeholder for now
#' }
get_snp_tree <- function(cellid_bam_table,
                         temp_dir = tempdir(),
                         slurm_base = paste(getwd(), "/slurmOut", sep = ""),
                         account = "gdrobertslab",
                         ploidy,
                         ref_fasta,
                         min_depth = 5,
                         submit = TRUE) {
    check_cellid_bam_table(cellid_bam_table)

    # Check that conda environment py3_10 exists, and if not, create it
    ## Check that the conda command is available
    check_cmd("conda")
    if (system("conda env list | grep py3_10") != 0) {
        message("Creating required conda environment py3_10")
        system("conda env create -f environment.yml py3_10.yml")
    }

    # Warn that this is going to take a while
    message("Hold onto your hat and get a coffee, this will take a while.")

    bam_files <- unique(cellid_bam_table$bam_file)
    # Write out the table if cell ids, groups and bam files
    parallel::mclapply(bam_files,
                       function(x) {
                            sub_cellid_bam_table <-
                                cellid_bam_table %>%
                                dplyr::filter(bam_file == x) %>%
                                dplyr::select(cell_barcode, cell_group)

                            bam_stub <- stringr::str_remove(x, ".bam")
                            call_snps(cellid_bam_table = sub_cellid_bam_table,
                                      bam_to_use = x,
                                      sam_dir = paste0(temp_dir,
                                                       "split_sams/",
                                                       bam_stub),
                                      bcf_dir = paste0(temp_dir,
                                                       "split_bcfs/",
                                                       bam_stub),
                                      slurm_base = slurm_base,
                                      account = account,
                                      ploidy = ploidy,
                                      ref_fasta = ref_fasta,
                                      submit = submit)
                       },
                       mc.cores = length(bam_files))

    # The output from the previous step is a folder for each bam file located
    # in temp_dir/split_bcfs/. Next merge all the bcf files
    merge_bcfs(bcf_folder = temp_dir,
               out_file = paste0(temp_dir, "/merged.bcf"),
               submit = submit)

    # Calculate a distance matrix
    dist_from_bcf(bcf = paste0(temp_dir, "/merged.bcf"),
                  out_file = paste0(temp_dir, "/distances.txt"),
                  submit = submit)

    # Read in the distance matrix and make a tree
    tree_out <- calc_tree(bcf = paste0(temp_dir, "/merged.bcf"))

    return(tree_out)
}

#' Check that the cellid_bam_table is of the proper format
#'
#' @param cellid_bam_table the cellid_bam_table to check
#'
#' @return 0 if the table is of the proper format, otherwise an error is thrown
check_cellid_bam_table <- function(cellid_bam_table) {
    # Columns should be cell_barcode, cell_group and bam_file
    if (any(!c("cell_barcode",
               "cell_group",
               "bam_file") %in% colnames(cellid_bam_table))) {
        stop("Columns should include cell_id, cell_group and bam_file")
    }

    # check that cell barcode is of format `[ATGC]+-1`
    if (!all(stringr::str_detect(cellid_bam_table$cell_barcode,
                                 "^[ATGC]+-1$"))) {
        warning("!!!\nCell barcodes should be of format `^[ATGC]+-1$`.\n",
                "We're going to try anyways, but this might not work.\n!!!")
    }

    # check that each cell_group is unique to a single bam_file
    # if the same group label is present in multiple bam_files, when we call
    # snps and merge the bcfs, we'll get duplicate column names and errors
    bams_per_group <-
        cellid_bam_table %>%
        dplyr::select(-cell_barcode) %>%
        dplyr::distinct() %>%
        dplyr::group_by(cell_group) %>%
        dplyr::summarize(n = dplyr::n()) %>%
        dplyr::pull(n) %>%
        max()
    if (bams_per_group != 1) {
        stop("Each cell_group should be unique to a single bam file")
    }
    return(0)
}

# ploidy options: GRCh37, GRCh38, X, Y, 1, mm10_hg19, mm10

call_snps <- function(cellid_bam_table,
                      bam_to_use,
                      sam_dir,
                      bcf_dir,
                      slurm_base = paste0(getwd(), "/slurmOut_call-%j.txt"),
                      account = "gdrobertslab",
                      ploidy,
                      ref_fasta,
                      min_depth = 5,
                      submit = TRUE,
                      cleanup = TRUE) {
    # Check that the sam and bcf directories exist, and if not, create them
    if (!dir.exists(sam_dir)) {
        dir.create(sam_dir)
    }
    if (!dir.exists(bcf_dir)) {
        dir.create(bcf_dir)
    }
    # write out the cell ids to a file with two columns: cell_id, cell_group
    cell_file <- paste0(sam_dir, "/cell_ids.txt")
    readr::write_tsv(dplyr::select(cellid_bam_table, cell_barcode, cell_group),
                     file = cell_file,
                     col_names = FALSE)

    # call getBarcodesFromBam.py on the bam file and the cell id file by reading
    # in a template and substituting out the placeholder fields
    # write the bam files to a subfolder for each source bam
    py_file <-
        paste0(find.package("rrrSingleCellUtils"),
               "/exec/getBarcodesFromBam.py")

    replace_tibble_split <-
        dplyr::tribble(
            ~find,                      ~replace,
            "placeholder_account",      account,
            "placeholder_slurm_out",    slurm_base,
            "placeholder_cell_file",    cell_file,
            "placeholder_bam_file",     bam_to_use,
            "placeholder_sam_dir",      paste0(sam_dir, "/"),
            "placeholder_py_file",      py_file
        )

    result <-
        use_sbatch_template(replace_tibble_split,
                            "snp_call_splitbams_template.job",
                            warning_label = "Bam splitting",
                            submit = submit)

    ploidy <- pick_ploidy(ploidy)

    array_max <-
        list.files(path = sam_dir, pattern = ".sam") %>%
        length()  - 1

    replace_tibble_snp <-
        dplyr::tribble(
            ~find,                      ~replace,
            "placeholder_account",      account,
            "placeholder_slurm_out",    slurm_base,
            "placeholder_array_max",    as.character(array_max),
            "placeholder_sam_dir",      sam_dir,
            "placeholder_ref_fasta",    ref_fasta,
            "placeholder_ploidy",       ploidy,
            "placeholder_bam_file",     bam_to_use,
            "placeholder_bcf_dir",      bcf_dir,
            "placeholder_min_depth",    as.character(min_depth)
        )

    # Call mpileup on each split_bams folder using a template and substituting
    # out the placeholder fields and index the individual bcf files
    result <-
        use_sbatch_template(replace_tibble_snp,
                            "snptree_call_snp_template.txt",
                            warning_label = "Calling SNPs",
                            submit = submit)

    # Delete contents of the split_sams folder
    if (cleanup) {
        unlink(sam_dir, recursive = TRUE)
    }
    return(0)
}

#' Transform the ploidy argument into a valid argument for bcftools
#'
#' @param ploidy Either a file path to a ploidy file, or a string indicating
#'  the ploidy.
#' @return A string that can be passed to use_sbatch_template() to fill in
#'  placeholder_ploidy
#' @details GRCh37 is hg19, GRCh38 is hg38, X, Y, 1, mm10_hg19 is our mixed
#' species reference with species prefixes on chromosomes, mm10 is mm10
pick_ploidy <- function(ploidy) {
    if (file.exists(ploidy)) {
        return(paste("--ploidy-file", ploidy))
    } else if (file.exists(paste0(find.package("rrrSingleCellUtils"),
                                  "/extdata/",
                                  ploidy,
                                  "_ploidy.txt"))) {
        return(paste0("--ploidy-file ",
                      find.package("rrrSingleCellUtils"),
                      "/extdata/",
                      ploidy,
                      "_ploidy.txt"))
    } else if (ploidy %in% c("GRCh37", "GRCh38", "X", "Y", "1")) {
        return(paste("--ploidy", ploidy))
    } else {
        stop("Ploidy argument not valid. Did you mean to pass a file path or ",
             "spell something wrong?")
    }
}

merge_bcfs <- function(bcf_folder,
                       out_file) {
    # use template to merge bcfs and write out a distance matrix, substituting
    # out the placeholder fields

    # index the bcf file

    # remove individual bcf files
}

dist_from_bcf <- function(bcf,
                          out_file) {
    # use template to calculate a distance matrix, substituting out the
    # placeholder fields

    # write out the distance matrix
}

calc_tree <- function(bcf) {
    # read in the distance matrix and make a tree

    # Return a hclust tree
}
