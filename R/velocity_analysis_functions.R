#' Wrapper to Make Loom Files in R
#'
#' @param sobj Seurat object you want to run velocity analysis
#' @param loom_dir Folder to output loom files to. By default will create a
#' directory in your working directory called loom_files
#' @param id_col Name of metadata column marking what sample a given cell is
#' from.
#' @param species String communicating what species the cells are from.
#' Value should be either "human" or "mouse".
#' @param bam_paths A named list or vector containing the bam files for each
#' unique sample in id_col. Names should be the sample name and should match the
#' unique values in sobj[[id_col]]
#' @param cluster_account Your Franklin cluster user ID.
#' @param slurm_base The directory to write slurm output files to.
#' @param sbatch_base The prefix to use with the sbatch job file.
#'
#' @details This is a wrapper function around velocyto's "run" command. It will
#' take your seurat object and generate loom files using only the cells present
#' in your seurat object. The loom files can then be used to run velocity
#' analysis using scVelo.
#' 
#' @return A loom file for each unique ID present in id_col, output to loom_dir.
#' 
#' @export

r_make_loom_files <- function(sobj,
                              loom_dir = "loom_files/",
                              id_col = NULL,
                              species,
                              bam_paths,
                              cluster_account,
                              slurm_base = paste0(getwd(), "/slurmOut"),
                              sbatch_base = "sbatch_") {
    #make loom_dir end with a /
    loom_dir <- ifelse(endsWith(loom_dir, "/"),
                       loom_dir,
                       paste0(loom_dir, "/"))

    #get ids and store in a variable
    samp_ids <- unique(sobj[[id_col]])[, 1]

    #Make sure sample_id column exists
    if (is.null(id_col)) {
        stop(paste("No column with sample ID's provided.",
                   "Please specify what column these are found in."))
    }
    if (length(setdiff(samp_ids, names(bam_paths))) != 0) {
        stop(paste("Names of bam_paths and unique ids in sobj[[id_col]]",
                   "differ. Please ensure they match and try again. Exiting...")
        )
    }


    #get proper gtf file
    if (species == "human") {
        gtf_path <- "/home/gdrobertslab/lab/GenRef/10x-hg38/genes/genes.gtf.gz"
    } else if (species == "mouse") {
        gtf_path <- "/home/gdrobertslab/lab/GenRef/10x-mm10/genes/genes.gtf.gz"
    } else {
        stop(paste("Unknown species.",
                   "Species parameter must be either \"human\" or \"mouse\".",
                   "Exiting..."))
    }

    #Make local copy of gtf file and unzip
    if (!file.exists(paste0(species, "_genes.gtf"))) {
        system(paste0("cp ",
                      gtf_path,
                      " ",
                      species,
                      "_genes.gtf.gz; gunzip ",
                      species,
                      "_genes.gtf.gz"))
    }

    #Make directories for sbatch files and slurm output
    system("mkdir sbatch; mkdir sbatch/jobs; mkdir sbatch/output")
    
    #optionally create tmp_bcs
    if(!dir.exists("tmp_bcs")) dir.create("tmp_bcs")

    #optionally create tmp_bams
    if (!dir.exists("tmp_bams")) dir.create("tmp_bams")

    #Loop through ID's
    for (id in unique(samp_ids)) {
        #make temporary directory with barcodes for current sample
        bcs <- colnames(sobj)[sobj@meta.data[[id_col]] == id]
        #Remove any additional things added to barcode
        bcs <- unlist(str_extract_all(bcs, "[A T G C]{16}+\\-+[1-9]"))
        bc_path <- paste0("tmp_bcs/", id, ".tsv")
        write.table(bcs,
                    bc_path,
                    row.names = FALSE,
                    col.names = FALSE,
                    quote = FALSE)

        #get bam path
        bam_path <- bam_paths[[id]]
        #copy bam to local location
        tmp_bam_path <- paste0("tmp_bams/", id, ".bam")
        system(paste("cp", bam_path, tmp_bam_path))

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

        replace_tbl <-
            tribble(~find, ~replace,
                    "placeholder_account", cluster_account,
                    "placeholder_slurm_out", paste0("sbatch/output/", id),
                    "placeholder_cell_file", bc_path,
                    "placeholder_bam_file", tmp_bam_path,
                    "placeholder_loom_dir", loom_dir,
                    "placeholder_env_path", env_path,
                    "placeholder_gtf_file", paste0(species, "_genes.gtf"),
                    "placeholder_sampleid", id)

        use_sbatch_template(replace_tibble = replace_tbl,
                            template = paste0(rrrscu,
                                              "/make_loom_files.sh"),
                            submit = TRUE,
                            file_dir = "sbatch/jobs")
    }
    #Remove tmp_bcs, tmp_bams, and genes files
    system("rm -r tmp_bcs; rm -r tmp_bams; rm -f *_genes.gtf*")
}

#' Save Off Metadata for Velocity Analysis
#'
#' @param sobj Seurat object; must have a column titled sample_id
#' @param sobj_id The specific seurat object ID for which you want metadata;
#' required even if only one sample in sobj
#' @param output_dir The directory you wish to save your metadata to
#' @param vars_to_keep Metadata columns you want saved off along with sample_id
#' and UMAP and PCA embeddings
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
                         vars_to_keep = NULL) {

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
        tmp_md <- dplyr::select(tmp_ob@meta.data, any_of(vars_to_keep)) %>%
            rownames_to_column("bc")

        #save off reduction coordinates if present
        for (red in names(tmp_ob@reductions)) {
            tmp_md <- cbind(tmp_md,
                            Seurat::Embeddings(tmp_ob, reduction = red))
        }

        #change rownames so they match format in the loom files
        tmp_md$bc <-
            paste0(id,
                   ":",
                   unlist(stringr::str_extract_all(tmp_md$bc, "[A T G C]{16}")),
                   "x")

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