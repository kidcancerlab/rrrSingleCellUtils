#' Generate a tree using SNPs derived from single cell clusters
#'
#' @param cellid_bam_table A tibble with three columns: cell_barcode,
#'  cell_group and bam_file. The cell_barcode column should contain the cell
#'  barcode, the cell_group column should contain the cluster label and the
#'  bam_file column should contain the path to the bam file for that cell.
#' @param temp_dir The directory to write temporary files to.
#' @param slurm_base The directory to write slurm output files to.
#' @param sbatch_base The prefix to use with the sbatch job file.
#' @param account The hpc account to use.
#' @param ploidy Either a file path to a ploidy file, or a string indicating
#'  the ploidy.
#' @param ref_fasta The path to the reference fasta file.
#' @param min_depth The minimum depth to use when calling SNPs.
#' @param min_sites_covered The minimum number of sites that must be covered by
#'  a cell_group to be included in the tree.
#' @param submit Whether to submit the sbatch jobs to the cluster or not.
#'
#' @details for ploidy, GRCh37 is hg19, GRCh38 is hg38, X, Y, 1, mm10_hg19 is
#'  our mixed species reference with species prefixes on chromosomes, mm10 is
#'  mm10
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
                         slurm_base = paste0(getwd(), "/slurmOut"),
                         sbatch_base = "sbatch_",
                         account = "gdrobertslab",
                         ploidy,
                         ref_fasta,
                         min_depth = 5,
                         min_sites_covered = 500,
                         submit = TRUE,
                         cleanup = TRUE) {
    check_cellid_bam_table(cellid_bam_table)

    ## Check that the conda command is available
    check_cmd("conda")
    # Check that conda environment rrrSingleCellUtils_py3_10 exists, and if not,
    # create it
    confirm_conda_env()

    # Check that temp_dir exists, and if not, create it
    if (!dir.exists(temp_dir)) {
        message("Creating temporary directory: ", temp_dir)
        dir.create(temp_dir)
    }

    # Warn that this is going to take a while
    message("Hold onto your hat and get a coffee, this will take a while.")

    bam_files <- unique(cellid_bam_table$bam_file)

    results <-
        parallel::mclapply(seq_len(length(bam_files)),
                           function(i) {
            bam_name <- bam_files[i]

            sub_cellid_bam_table <-
                cellid_bam_table %>%
                dplyr::filter(bam_file == bam_name) %>%
                dplyr::select(cell_barcode, cell_group)

            call_snps(cellid_bam_table = sub_cellid_bam_table,
                      bam_to_use = bam_name,
                      sam_dir = paste0(temp_dir,
                                       "/split_sams_",
                                       i,
                                       "/"),
                      bcf_dir = paste0(temp_dir,
                                       "/split_bcfs_",
                                       i,
                                       "/"),
                      slurm_base = paste0(slurm_base, "_mpileup-%j.out"),
                      sbatch_base = sbatch_base,
                      account = account,
                      ploidy = ploidy,
                      min_depth = min_depth,
                      ref_fasta = ref_fasta,
                      submit = submit,
                      cleanup = cleanup)
            },
                        mc.cores = length(bam_files))

    # The output from the previous step is a folder for each bam file located
    # in temp_dir/split_bcfs/. Next merge all the bcf files and calculate a
    # distance matrix using my slow python script

    merge_bcfs(bcf_in_dir = paste0(temp_dir,
                                   "/split_bcfs_[0-9]*/"),
               out_bcf = paste0(temp_dir, "/merged.bcf"),
               out_dist = paste0(temp_dir, "/distances"),
               submit = submit,
               slurm_out = paste0(slurm_base, "_merge-%j.out"),
               sbatch_base = sbatch_base,
               account = account,
               cleanup = cleanup)

    # Read in the distance matrix and make a tree
    dist_tree <-
        calc_tree(matrix_file = paste0(temp_dir, "/distances.tsv"),
                  counts_file = paste0(temp_dir, "/distances_tot_count.tsv"),
                  min_sites = min_sites_covered)

    return(dist_tree)
}

#' Use the output from get_snp_tree to label cells
#'
#' @param sobject A Seurat object
#' @param dist_tree A hclust tree output from get_snp_tree()
#' @param group_col_name The name of the column in the Seurat object that was
#'  used when assigning cells to groups for get_snp_tree()
#' @param control_group The name of the control group. If this is not NULL,
#'  then the group containing the control will be labeled "normal" and the
#'  other group(s) will be labeled "tumor".
#' @param cut_n_groups The number of groups to cut the tree into.
#' @param cut_dist The distance to cut the tree at. This can be derived from the
#'  y-axis of plotting the output from get_snp_tree().
#'
#' @return A Seurat object with a new columns in the metadata slot called
#'  tree_group and, if control_group is not NULL, snp_tumor_call.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'    placeholder here for now
#' }
label_tree_groups <- function(sobject,
                              dist_tree,
                              group_col_name = "used_clusters",
                              normal_groups,
                              cut_n_groups = NULL,
                              cut_dist = NULL) {
    tree_groups <-
        stats::cutree(dist_tree,
                      k = cut_n_groups,
                      h = cut_dist) %>%
        tibble::enframe(name = group_col_name,
                        value = "tree_group") %>%
        dplyr::mutate(tree_group = LETTERS[tree_group])

    # if tree_group contains the control, then rename it to "normal" and label
    # the other group(s) as "tumor"
    if (!is.null(normal_groups)) {
        tree_groups <-
            tree_groups %>%
            dplyr::group_by(tree_group) %>%
            dplyr::mutate(snp_tumor_call = ifelse(any(normal_groups %in%
                                                   get(group_col_name)),
                                                  "normal",
                                                  "tumor"))
    }
    sobj_out <-
        sobject@meta.data %>%
        tibble::rownames_to_column("cell_barcode") %>%
        dplyr::left_join(tree_groups, by = group_col_name) %>%
        dplyr::select(dplyr::all_of(c(colnames(tree_groups),
                                      "cell_barcode"))) %>%
        tibble::column_to_rownames("cell_barcode") %>%
        Seurat::AddMetaData(object = sobject)

    return(sobj_out)
}

#' Confirm conda environment and create it if it doesn't exist
confirm_conda_env <- function() {
    if (system("conda env list | grep rrrSingleCellUtils_py3_10",
               ignore.stdout = TRUE) != 0) {
        conda_yml_file <-
            paste0(find.package("rrrSingleCellUtils"),
                   "/py3_10.yml")
        message("Creating required conda environment rrrSingleCellUtils_py3_10")
        system(paste0("conda env create -n rrrSingleCellUtils_py3_10 --file ",
                      conda_yml_file))
    }
    return()
}

#' Check that the cellid_bam_table input is of the proper format
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

    # Check that there is data in the table
    if (nrow(cellid_bam_table) == 0) {
        stop("No data in cellid_bam_table")
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
        dplyr::summarize(n_bams = dplyr::n()) %>%
        dplyr::pull(n_bams) %>%
        max()
    if (bams_per_group != 1) {
        stop("Each cell_group should be unique to a single bam file")
    }

    # check that there are no spaces in either the cell_group or bam_file
    if (any(stringr::str_detect(cellid_bam_table$cell_group, " ")) ||
        any(stringr::str_detect(cellid_bam_table$bam_file, " "))) {
        stop("Neither cell_group nor bam_file should not contain spaces")
    }
    return(0)
}

#' Call SNPs for a single bam file
#'
#' @param cellid_bam_table A table with columns cell_id, cell_group and bam_file
#' @param bam_to_use The bam file to use
#' @param sam_dir The directory to write the sam files to
#' @param bcf_dir The directory to write the bcf files to
#' @param slurm_base The base name for the slurm output files. Don't include path
#' @param account The hpc account to use in slurm scripts
#' @param ploidy The ploidy of the organism. See details for more information
#' @param ref_fasta The reference fasta file to use
#' @param min_depth The minimum depth to use when calling snps
#' @param submit Whether or not to submit the slurm scripts
#' @param cleanup Whether or not to clean up the sam files afterwards
#'
#' @return 0 if the snps were called successfully
#'
#' @details GRCh37 is hg19, GRCh38 is hg38, X, Y, 1, mm10_hg19 is our mixed
#'
call_snps <- function(cellid_bam_table,
                      bam_to_use,
                      sam_dir,
                      bcf_dir,
                      slurm_base = paste0(getwd(), "/slurmOut_call-%j.txt"),
                      sbatch_base = "sbatch_",
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
                     col_names = FALSE,
                     progress = FALSE)

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
                            submit = submit,
                            file_dir = ".",
                            temp_prefix = paste0(sbatch_base, "split_"))

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
                            "snp_call_mpileup_template.job",
                            warning_label = "Calling SNPs",
                            submit = submit,
                            file_dir = ".",
                            temp_prefix = paste0(sbatch_base, "mpileup_"))

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
        warning("Ploidy argument not valid. Did you mean to pass a file path or ",
                "spell something wrong?",
                immediate. = TRUE)
        stop()
    }
}

#' Merge bcfs generated by call_snps() and write out a distance matrix
#'
#' @param bcf_in_dir The directory containing the bcfs to merge
#' @param out_bcf The path to the output bcf file
#' @param out_dist The path to the output distance matrix
#' @param submit Whether to submit the job to slurm
#' @param account The cluster account to use in the slurm script
#' @param slurm_out The name of the slurm output file
#' @param cleanup Whether to delete the bcfs after merging
#'
#' @return 0 if successful
merge_bcfs <- function(bcf_in_dir,
                       out_bcf,
                       out_dist,
                       submit = TRUE,
                       account = "gdrobertslab",
                       slurm_out = "slurmOut_merge-%j.txt",
                       sbatch_base = "sbatch_",
                       cleanup = TRUE) {
    py_file <-
        paste0(find.package("rrrSingleCellUtils"),
               "/exec/vcfToMatrix.py")

    # use template to merge bcfs and write out a distance matrix, substituting
    # out the placeholder fields
    replace_tibble_merge_dist <-
        dplyr::tribble(
            ~find,                      ~replace,
            "placeholder_account",      account,
            "placeholder_slurm_out",    slurm_out,
            "placeholder_bcf_out",      out_bcf,
            "placeholder_py_file",      py_file,
            "placeholder_bcf_dir",      bcf_in_dir,
            "placeholder_dist_out",     out_dist
        )

    # Call mpileup on each split_bams folder using a template and substituting
    # out the placeholder fields and index the individual bcf files
    result <-
        use_sbatch_template(replace_tibble_merge_dist,
                            "snp_call_merge_dist_template.job",
                            warning_label = "Calling SNPs",
                            submit = submit,
                            file_dir = ".",
                            temp_prefix = paste0(sbatch_base, "merge_"))

    # remove individual bcf files
    if (cleanup) {
        unlink(bcf_in_dir, recursive = TRUE)
    }
    return(0)
}

#' Calculate a tree from a distance matrix output from merge_bcfs()
#'
#' @param matrix_file The path to the distance matrix output from merge_bcfs().
#' @param counts_file The path to the counts file output from merge_bcfs().
#' @param min_sites The minimum number of sites that must be covered by a
#'  cell_group to be included in the tree.
#'
#' @export
#'
#' @return A hclust object
calc_tree <- function(matrix_file,
                      counts_file,
                      min_sites = 500) {
    # filter out any clusters with too few SNP sites
    high_counts <-
        read.table(counts_file,
                   header = TRUE,
                   row.names = 1,
                   sep = "\t") %>%
        as.matrix() %>%
        diag() %>%
        purrr::keep(~.x > min_sites) %>%
        names()

    # read in the distance matrix and make a tree
    dist_tree <-
        read.table(matrix_file,
                   header = TRUE,
                   row.names = 1,
                   sep = "\t") %>%
        tibble::rownames_to_column("sample_1") %>%
        tidyr::pivot_longer(names_to = "sample_2",
                            values_to = "snp_dist",
                            -sample_1) %>%
        dplyr::filter(!is.na(snp_dist) &
                      sample_1 %in% high_counts &
                      sample_2 %in% high_counts) %>%
        dplyr::group_by(sample_1) %>%
        dplyr::mutate(group_count = dplyr::n()) %>%
        dplyr::filter(group_count > 1) %>%
        dplyr::select(-group_count) %>%
        tidyr::pivot_wider(names_from = sample_2,
                           values_from = snp_dist) %>%
        tibble::column_to_rownames("sample_1") %>%
        as.matrix() %>%
        stats::dist() %>%
        stats::hclust()

    # Return a hclust tree
    return(dist_tree)
}

#' Use the metadata in a Seurat object to determine which clusters are
#' predominantly control cell types
#'
#' @param sobject A Seurat object
#' @param normal_celltypes A vector of cell types that are considered non-tumor
#' @param cluster_col The name of the column that was used to define clusters in
#'  get_snp_tree()
#' @param celltype_col The name of the column that has cell type information
#' @param min_prop_control The minimum proportion of cells in a cluster that
#'  are control cell types to consider that cluster a control cluster
#'
#' @return A vector of cluster names that are predominantly control cell types
#'
#' @export
match_celltype_clusters <- function(sobject,
                                    normal_celltypes,
                                    cluster_col,
                                    celltype_col,
                                    min_prop_control = 0.5) {

    control_clusters <-
        sobject@meta.data %>%
        dplyr::select(dplyr::all_of(c(cluster_col, celltype_col))) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(cluster_col))) %>%
        dplyr::filter(sum(get(celltype_col) %in%
                                                   control_celltypes) /
                                               dplyr::n() > min_prop_control) %>%
        dplyr::pull(cluster_col) %>%
        unique()

    return(control_clusters)
}
