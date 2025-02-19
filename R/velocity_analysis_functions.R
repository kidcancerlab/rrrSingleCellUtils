#' Wrapper to Make Loom Files in R
#'
#' @param input_table Dataframe or tibble with columns titled "sample_id",
#' "h5_path", "bam_path", and "gtf_path"
#' @param out_dir Folder to output loom files and sbatch output to. By default
#' will create a directory in your working directory called loom_output/samples.
#' Each sample_id found in input_table will have its own directory within
#' out_dir
#' @param cluster_account Your Franklin cluster user ID.
#' @param slurm_base The directory to write slurm output files to.
#' @param sbatch_base The prefix to use with the sbatch job file.
#'
#' @details This is a wrapper function around velocyto's "run" command. It will
#' take a table of sid's and other info and create loom files based on all cells
#' that are in the h5 object. The loom files can then be used to run velocity
#' analysis using scVelo. Slurm error and std output will be written to
#' out_dir/sample_id/.
#'
#' @return A loom file for each unique ID present in id_col, output to loom_dir.
#'
#' @export


r_make_loom_files <- function(input_table,
                              out_dir = "loom_output/samples",
                              cluster_account,
                              slurm_base = paste0(getwd(), "/slurmOut"),
                              sbatch_base = "sbatch_") {
    
    #Check for column names
    if (!all((c("sample_id",
                "h5_path",
                "bam_path",
                "gtf_path") %in% colnames(input_table)))) {
        missing_cols <- setdiff(c("sample_id", "h5_path", "bam_path", "gtf_path"), #nolint
                                colnames(input_table))
        stop(paste("Column(s)",
                   paste(missing_cols,
                         collapse = ", ")),
                   " not found in input table. Exiting...")
    }

    #Check all h5, bam, and gtf files exist
    if (!all(file.exists(input_table$h5_path))) {
        stop(paste("H5 files not found for",
                   input_table$sample_id[which(!file.exists(input_table$h5_path))])) #nolint
    }
    if (!all(file.exists(input_table$bam_path))) {
        stop(paste("bam files not found for",
                   input_table$sample_id[which(!file.exists(input_table$bam_path))])) #nolint
    }
    if (!all(file.exists(input_table$gtf_path))) {
        stop(paste("gtf files not found for",
                   input_table$sample_id[which(!file.exists(input_table$gtf_path))])) #nolint
    }

    rownames(input_table) <- input_table$sample_id

    #make out_dir end with a /
    out_dir <- ifelse(endsWith(out_dir, "/"),
                      out_dir,
                      paste0(out_dir, "/"))

    system(paste0("mkdir -p ", out_dir))

    for (sid in rownames(input_table)) {
        #make sample output folders
        sid_out <- paste0(out_dir, sid)
        system(paste("mkdir -p", sid_out))

        #Make local copy of gtf file and unzip
        if (endsWith(input_table[sid, ]$gtf_path, ".gz")) {
            new_gene_path <- paste0(sid_out, "/genes.gtf.gz")
            system(paste0("cp ",
                         input_table[sid, ]$gtf_path,
                         " ",
                         new_gene_path,
                         "; gunzip ",
                         new_gene_path,
                         " ", 
                         sid_out,
                         "/genes.gtf"
                         ))
        } else {
            system(paste0("cp ",
                         input_table[sid, ]$gtf_path,
                         " ",
                         sid_out,
                         "/genes.gtf"))
        }

        #read in h5 object
        h5_object <- Read10X_h5(input_table[sid, ]$h5_path)
        #make sure only get gene expression data
        if (class(h5_object) == "list") {
            h5_object <- h5_object[["Gene Expression"]]
        }

        #make temporary directory with barcodes for current sample
        bcs <- colnames(h5_object)
        write.table(bcs,
                    paste0(sid_out, "/tmp_bcs.tsv"),
                    row.names = FALSE,
                    col.names = FALSE,
                    quote = FALSE)

        #copy over bams
        old_bam_path <- input_table[sid, ]$bam_path
        tmp_bam_path <- paste0(sid_out, "/tmp_bam.bam")
        if (!file.exists(tmp_bam_path)) {
            system(paste("cp", old_bam_path, tmp_bam_path))
        }
    }
    #Make conda environment
    #going to make environment in location of R installation, it's likely
    #that people won't have any conda environments here
    #get location of rrrSingleCellUtils and append env name to it
    rrrscu <- find.package("rrrSingleCellUtils")
    env_path <- paste0(rrrscu, "/r_rna_velo")

    #check conda environment doesn't exist before creating
    conda_envs <- system("conda info --envs", intern = TRUE)

    if (sum(grepl(pattern = env_path, x = conda_envs)) == 0) {
        #only make conda environment if it doesn't already exist
        exists_conda <-
            system(paste0("conda env create -p ",
                          env_path,
                          " -f ",
                          paste0(rrrscu,
                                 "/make_environment.yml")))
    }

    #create bash array of sample_ids
    ids <- rownames(input_table)
    id_array <- paste(ids, collapse = " ")

    replace_tbl <-
        tibble::tribble(
            ~find,                      ~replace,
            "placeholder_account",      cluster_account,
            "placeholder_slurm_out",    paste0(out_dir, "sbatch/out_"),
            "placeholder_slurm_error",  paste0(out_dir, "sbatch/error_"),
            "placeholder_env_path",     env_path,
            "placeholder_max_array",    as.character(length(ids) - 1),
            "placeholder_id_array",     id_array,
            "placeholder_out_dir",      out_dir
        )

    #make temp dir for temp sbatch job files
    system(paste0("mkdir ", out_dir, "sbatch"))
    use_sbatch_template(replace_tibble = replace_tbl,
                        template = paste0(rrrscu, "/make_loom_files.sh"),
                        submit = TRUE,
                        file_dir = paste0(out_dir, "sbatch/"))
}

#' Save Off Metadata for Velocity Analysis
#'
#' @param sobj Seurat object; must have a column titled sample_id
#' @param sobj_id The specific seurat object ID for which you want metadata;
#' required even if only one sample in sobj
#' @param output_dir The directory you wish to save your metadata to
#' @param vars_to_keep Metadata columns you want saved off along with sample_id
#' and UMAP and PCA embeddings
#' @param handle_n_of_1 boolean saying whether or not you want to handle the
#' special case when a sample has only one observation
#'
#' @details This is a helper function for running RNA velocity analysis. While
#' your Seurat object may contain multiple samples, the loom files are
#' constructed for each sample_id present in that object, and the anndata
#' objects created from those loom files are for each individual sample_id.
#' The function assumes that you made the loom files using the shell script
#' "make_loom_files.sh", which changes the row names and appends a unique ID
#' to the end of each sample when saving it off.
#'
#' @return A csv for each sample ID containing the sample ID, the reductions,
#' and any user specified columns of metadata in output_dir
#'
#' @export

write_off_md <- function(sobj,
                         id_col,
                         output_dir,
                         vars_to_keep = NULL,
                         handle_n_of_1 = TRUE) {

    output_dir <- ifelse(endsWith(output_dir, "/"),
                         substr(output_dir, 1, nchar(output_dir) - 1),
                         output_dir)

    #Make sure sample_id column exists
    if (is.null(id_col)) {
        stop(paste("No column with sample ID's provided.",
                   "Please specify what column these are found in."))
    }

    #get ids and store in a variable
    samp_ids <- unique(sobj[[id_col]])[, 1]

    #add id_col to vars_to_keep
    vars_to_keep <- c(vars_to_keep, id_col)

    for (id in samp_ids) {
        #subset object for current id
        tmp_ob <- sobj[, sobj@meta.data[[id_col]] == id]
        tmp_md <-
            dplyr::select(
                tmp_ob@meta.data,
                dplyr::any_of(vars_to_keep)
            ) %>%
            tibble::rownames_to_column("bc")

        #save off reduction coordinates if present
        for (red in names(tmp_ob@reductions)) {
            tmp_md <- cbind(tmp_md,
                            Seurat::Embeddings(tmp_ob, reduction = red))
        }

        #change rownames so they match format in the loom files
        #account for case of n = 1
        if (nrow(tmp_md) == 1 && handle_n_of_1) {
            tmp_md$bc <-
                paste0(
                    id,
                    ":",
                    unlist(
                        stringr::str_extract_all(
                            tmp_md$bc,
                            "[ATGC]{16}\\-[0-9]"
                        )
                    )
                )
        } else {
            tmp_md$bc <-
                paste0(
                    id,
                    ":",
                    unlist(
                        stringr::str_extract_all(
                            tmp_md$bc,
                            "[A T G C]{16}"
                        )
                    ),
                    "x"
                )
        }

        #write off metadata file
        if (is.null(output_dir)) {
            utils::write.csv(tmp_md,
                      paste0(id, "_metadata.csv"),
                      row.names = FALSE)
        } else {
            if (!dir.exists(output_dir)) dir.create(output_dir)
            utils::write.csv(tmp_md,
                    paste0(output_dir, "/", id, "_metadata.csv"),
                    row.names = FALSE)
        }
    }
}

#' Use a Sbatch Template to Submit a Job to the Cluster
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
        readr::read_file(template)

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
        return_val <- system(paste("sbatch", temp_file))
    } else {
        return_val <- 0
    }

    if (return_val != 0) {
        stop(paste0(warning_label,
                    " sbatch submission failed. Error code ",
                    return_val))
    }
    return(0)
}
