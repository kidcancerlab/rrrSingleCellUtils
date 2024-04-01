utils::globalVariables(c(
    ".",
    "Abort",
    "avg.exp.scaled",
    "avg_log2FC",
    "bam_file",
    "CB",
    "cell_barcode",
    "cell_group",
    "cid",
    "cluster",
    "Cluster",
    "col.fill",
    "control_celltypes",
    "database",
    "download_data",
    "exp_type",
    "fastq_folder_suffix",
    "fastqs",
    "feature",
    "features.plot",
    "freq",
    "Freq",
    "from",
    "gene",
    "group_count",
    "id",
    "junk",
    "label",
    "label14",
    "label30",
    "library_type",
    "ligand_target_matrix",
    "ligands",
    "ligands_bona_fide",
    "link_folder",
    "lr_network",
    "lr_network_strict",
    "lt_14",
    "lt_30",
    "max_val",
    "median_val",
    "method",
    "min_val",
    "n_bams",
    "num_clusters",
    "p_val_adj",
    "pearson",
    "pct.exp",
    "Phase",
    "Proportion",
    "Protocol",
    "rainbow",
    "receptors",
    "receptors_bona_fide",
    "res_vals",
    "rm_path",
    "run_cellranger_count",
    "run_cellranger_mkfastq",
    "sample_1",
    "sample_2",
    "Sample_ID",
    "Sample_Project",
    "score",
    "sd_val",
    "sil_vals",
    "sil_width",
    "snp_dist",
    "suffix",
    "sobj",
    "tar_folder",
    "test_ligand",
    "to",
    "to_make_sobj",
    "tree_group",
    "tx_id",
    "tx_name",
    "value",
    "Var1",
    "weight",
    "weighted_networks",
    "weighted_networks_lr",
    "x",
    "y"
    ))


#' Use a sbatch template to submit a job to the cluster
#'
#' @param replace_tibble A tibble with two columns, find and replace.
#'  The find column should contain the placeholder text to be replaced and the
#'  replace column should contain the text to replace it with.
#' @param template The name of the template file to use. This should be a file
#'  in the rrrSingleCellUtils/inst folder.
#' @param file_dir The directory to write the temporary sbatch file to.
#' @param temp_ext The extension to use for the temporary sbatch file.
#' @param temp_prefix The prefix to use for the temporary sbatch file.
#' @param warning_label A string to use in the warning message if the sbatch
#'  submission fails.
#' @param submit Whether to actually submit the sbatch job or just write the
#' sbatch file. If FALSE, the sbatch file will be written but not submitted.
#'
#' @return 0 if the sbatch submission was successful, otherwise an error is
#'  thrown
use_sbatch_template <- function(replace_tibble,
                                template,
                                file_dir = tempdir(),
                                temp_ext = ".sh",
                                temp_prefix = "sbatch_",
                                warning_label = "",
                                submit = TRUE) {
    sbatch_template <-
        readr::read_file(paste0(find.package("rrrSingleCellUtils"),
                                "/",
                                template))

    # Replace placeholders with real data
    for (i in seq_len(nrow(replace_tibble))) {
        sbatch_template <-
        stringr::str_replace_all(sbatch_template,
                                 pattern = replace_tibble$find[i],
                                 replacement = replace_tibble$replace[i])
    }

    temp_file <-
        tempfile(fileext = temp_ext,
                 tmpdir = file_dir,
                 pattern = temp_prefix)

    readr::write_file(sbatch_template, file = temp_file)

    if (submit == TRUE) {
        return_value <- system(paste("sbatch", temp_file))
    } else {
        return_value <- 0
    }

    if (return_value != 0) {
        big_problem(paste0(warning_label,
                           " sbatch submission failed. Error code ",
                           return_value,
                           ". Check the sbatch file ",
                           temp_file))
    }

    if (return_value == 0) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Check that a command is available on the system
#'
#' @param cmd The command to check for
#' @return 0 if the command is available, otherwise an error is thrown
check_cmd <- function(cmd) {
    if (Sys.which(cmd) == "") {
        stop(paste(cmd,
                   "command not found, do you need to load a module or add",
                   " this to your PATH?"))
    }
    return(0)
}