two_sobj <- qs::qread("/home/gdrobertslab/mvc002/analyses/roberts/dev/testSnps/output/patient_met_1/two_sobj.qs")

c_b_t <-
    two_sobj@meta.data %>%
    select(used_clusters) %>%
    rename(cell_group = used_clusters) %>%
    rownames_to_column("cell_barcode") %>%
    as_tibble() %>%
    mutate(bam_file = paste0("/home/gdrobertslab/lab/Counts/",
                        str_remove(cell_barcode, "_.*"),
                        "/outs/possorted_genome_bam.bam"),
            cell_barcode = str_remove(cell_barcode, ".+_"))

get_snp_tree <- function(cellid_bam_table,
                         temp_dir = tempdir(),
                         slurm_base = paste(getwd(), "/slurmOut", sep = ""),
                         account = "gdrobertlab",
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
                            bam_stub <- stringr::str_remove(x, ".bam")
                            call_snps(cellid_bam_table = cellid_bam_table,
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

check_cellid_bam_table <- function(cellid_bam_table) {
    # Columns should be cell_barcode, cell_group and bam_file
    if (any(!colnames(cellid_bam_table) %in% c("cell_barcode",
                                               "cell_group",
                                               "bam_file"))) {
        stop("Columns should be cell_id, cell_group and bam_file")
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
        dplyr::summarize(n = n()) %>%
        dplyr::pull(n) %>%
        max()
    if (bams_per_group != 1) {
        stop("Each cell_group should be unique to a single bam file")
    }
}

call_snps <- function(cellid_bam_table,
                      bam_to_use,
                      sam_dir,
                      bcf_dir,
                      slurm_base = paste(getwd(), "/slurmOut_call", sep = ""),
                      account = "gdrobertlab",
                      ploidy,
                      ref_fasta,
                      min_depth = 5,
                      submit = TRUE) {
    # filter the table using bam_to_use to a single bam file and
    # write out the cell ids to a file with two columns: cell_id, cell_group
    cell_file <- paste0(sam_dir, "/cell_ids.txt")
    cellid_bam_table %>%
        filter(bam_file == bam_to_use) %>%
        select(cell_barcode, cell_group) %>%
        write_tsv(cell_file,
                  col_names = FALSE)

    # call getBarcodesFromBam.py on the bam file and the cell id file by reading
    # in a template and substituting out the placeholder fields
    # write the bam files to a subfolder for each source bam
    py_file <-
        paste0(find.package("rrrSingleCellUtils"),
               "/exec/getBarcodesFromBam.py")

    replace_tibble_split <-
        tribble(~find,                      ~replace,
                "placeholder_account",      account,
                "placeholder_slurm_out",    slurm_base,
                "placeholder_cell_file",    cell_file,
                "placeholder_bam_file",     bam_to_use,
                "placeholder_sam_dir",      sam_dir,
                "placeholder_py_file",      py_file)

    result <-
        use_sbatch_template(replace_tibble_split,
                            "snp_call_splitbams_template.job",
                            warning_label = "Bam splitting",
                            submit = submit)

    ploidy <- pick_ploidy(ploidy)

    replace_tibble_snp <-
        tribble(~find,                      ~replace,
                "placeholder_account",      account,
                "placeholder_slurm_out",    slurm_base,
                "placeholder_array_max",    array_max,
                "placeholder_sam_dir",      sam_dir,
                "placeholder_ref_fasta",    ref_fasta,
                "placeholder_ploidy",       ploidy,
                "placeholder_bam_file",     bam_to_use,
                "placeholder_bcf_dir",      bcf_dir,
                "placeholder_min_depth",    min_depth)

    # Call mpileup on each split_bams folder using a template and substituting
    # out the placeholder fields and index the individual bcf files
    result <-
        use_sbatch_template(replace_tibble,
                            "snptree_call_snp_template.txt",
                            warning_label = "Calling SNPs",
                            submit = submit)

    # Delete contents of the split_sams folder
}

#' Transform the ploidy argument into a valid argument for bcftools
#'
#' @param ploidy Either a file path to a ploidy file, or a string indicating
#'  the ploidy.
#' @return A string that can be passed to use_sbatch_template() to fill in
#'  placeholder_ploidy
pick_ploidy <- function(ploidy) {
    if (file.exists(ploidy)) {
        return(paste("--ploidy-file", ploidy))
    } else if (file.exists(paste0(find.package("rrrSingleCellUtils"),
                                  "/extdata/",
                                  ploidy,
                                  ".txt"))) {
        return(paste0("--ploidy-file ",
                      find.package("rrrSingleCellUtils"),
                      "/extdata/",
                      ploidy,
                      ".txt"))
    } else if(ploidy %in% c("GRCh37", "GRCh38", "X", "Y", "1")) {
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
