#' Download and pre-process raw 10X data
#'
#' @param sample_info File containing sample info (see details)
#' @param domain Where the raw files are hosted
#' @param email Email for Slurm notifications
#' @param include_introns Should intronic reads be included in counts?
#' @param plots Should plots be generated?
#' @param slurm_base Base name for slurm output files
#' @param bcl_folder Path to write BCL files
#' @param fastq_folder Path to write fastq files
#' @param counts_folder Path to write counts files
#' @param ref_folder Path to 10x reference folders
#' @param sobj_folder Path to write Seurat objects
#' @param proc_threads Number of threads to use for processing
#'
#' @details This is a wrapper for several functions to get and process single
#'  cell data from the NCH IGM core. I have built the defaults to be specific to
#'  the Roberts lab, so you may need to carefully change the defaults if you
#'  want to use it outside of this context. The input is a sample sheet that
#'  outlines the information about each sample. The data are then downloaded
#'  using smbclient and then md5sum checked and untar’d.
#'  The data are then processed with cellranger mkfastq and either cellranger
#'  count or cellranger-dna cnv (depending on the exp_type argument). A seurat
#'  object is generated as well and filtered using auto_subset().
#'
#' @export
#'
#' @examples
#' \dontrun{
#' process_raw_data(sample_info = "testSampleInfoSheet_test.txt",
#' email = "matthew.cannon@nationwidechildrens.org")
#' }
process_raw_data <- function(sample_info,
                             domain = "//igmdata/igm_roberts",
                             email = "",
                             include_introns = FALSE,
                             plots = FALSE,
                             slurm_base = paste(getwd(), "/slurmOut", sep = ""),
                             bcl_folder = "/home/gdrobertslab/lab/BCLs_2",
                             fastq_folder = "/home/gdrobertslab/lab/FASTQs_2",
                             counts_folder = "/home/gdrobertslab/lab/Counts_2",
                             ref_folder = "/home/gdrobertslab/lab/GenRef",
                             sobj_folder = "/home/gdrobertslab/lab/SeuratObj",
                             cutoff_hist_folder = "/home/gdrobertslab/lab/SeuratObj",
                             proc_threads = 10) {
    # Make sure sample_info has unix line endings
    system(paste0("dos2unix ", sample_info))

    # For each sample, check if files are present to show that data are already
    # 1. Processed into a saved Seurat object
    # 2. Processed through cellranger_count
    # 3. Processed through cellranger_mkfastq
    # 4. Downloaded and md5sum checked
    # Each level negates the need to do the following levels
    # I'm filtering out aborted samples and samples that aren't scRNA-seq_3prime
    # or multiomics for now since I don't have the code to handle others yet
    sample_data <-
        readr::read_delim(sample_info,
                          delim = "\t",
                          col_names = TRUE,
                          trim_ws = TRUE,
                          show_col_types = FALSE) %>%
        dplyr::filter(Abort == "N" &
                      exp_type %in% c("MATAC",
                                      "MGEX",
                                      "3GEX")) %>%
        add_data_status(sample_data,
                        bcl_folder = bcl_folder,
                        fastq_folder = fastq_folder,
                        counts_folder = counts_folder,
                        sobj_folder = sobj_folder)

    #################### Download the things that need downloaded
    to_download <-
        dplyr::filter(sample_data,
                      download_data == TRUE)

    # Sometimes we're going to have multiple download folders and multiple
    # experiment types (such as multiomics) that need to be handled separately
    # Therefore, we loop over each link_folder
    #temp_sample_info_sheets <- list()

    pw <- getPass::getPass("Password for smbclient: ")

    for (dl_link in unique(to_download$link_folder)) {
        # Download, untar and check md5 sums of raw data
        message("Getting raw data from ", dl_link, ".")

        ############################### Need to handle the case in which the link is dead

        tar_list <-
            get_raw_data(link_folder = dl_link,
                         dest_folder = dest_folder,
                         domain = domain,
                         .pw = pw)
    }

    ###
    # For now, I am assuming there is one tar file per link_folder, but this
    # needs to be changed to handle multiple tar files per link_folder later
    ###

    ############### mkfastq the things that need fastqs
    to_mkfastq <-
        dplyr::filter(sample_data,
                      run_cellranger_mkfastq == TRUE)

    # Get unique combination of tar_folder and exp_type to use with the
    # following loop
    tar_exps <-
        to_mkfastq %>%
        dplyr::select(link_folder, tar_folder, exp_type, Sample_Project) %>%
        dplyr::distinct()

    # Loop over each tar_folder and exp_type combination to run mkfastq
    # number of cores used by mclapply doesn't matter since it's just submitting
    # them to the cluster
    parallel::mclapply(seq_len(nrow(tar_exps)),
                       mc.cores = 100,
                       function(i) {
        tar_f <- tar_exps$tar_folder[i]
        # Write out subsampled sampleInfoSheet with data from one link_folder/exp_type combo
        temp_sample_info_sheet <-
            tempfile(fileext = ".tsv")
        dplyr::left_join(tar_exps[i, ], sample_data) %>%
            readr::write_tsv(file = temp_sample_info_sheet)

        # Fix sample sheet for all folders with BCL files
        message("Fixing sample sheet in ", tar_f, ".")
        fix_sample_sheet(sample_sheet = paste0(bcl_folder, "/",
                                               tar_exps$Sample_Project[i], "/",
                                               tar_exps$tar_folder[i], "/",
                                               "/SampleSheet.csv"),
                         sample_info_sheet = temp_sample_info_sheet)

        # Run cellranger mkfastq
        message("Submitting slurm command to create fastq",
                "files using cellranger mkfastq.")
        cellranger_mkfastq(sample_info = temp_sample_info_sheet,
                           email = email,
                           tar_folders = tar_f,
                           bcl_folder = bcl_folder,
                           fastq_folder = fastq_folder,
                           slurm_out = paste0(slurm_base,
                                             "_mkfastq-%j.out"))
    })

    ###################### Count the things that need counting
    # Run cellranger count if scRNA-seq_3prime
    to_count <-
        dplyr::filter(sample_data,
                      run_cellranger_count == TRUE)

    # Need to feed cellranger_count subsets of the sample_data sheet based on
    # exp_type and Sample_Project
    to_count_tbl_list <-
        to_count %>%
        dplyr::group_by(Sample_Project, exp_type) %>%
        dplyr::group_split()

    # The number of cores used here doesn't matter, since it's just submitting
    # the jobs to the cluster
    parallel::mclapply(to_count_tbl_list,
                       mc.cores = 100,
                       function(x) {
        message("Submtting slurm command to run cellranger count.\n",
                "Slurm messages output to", paste0(slurm_base, "_count-%j.out"))

        cellranger_count(sample_info = x,
                         email = email,
                         counts_folder = counts_folder,
                         fastq_folder = fastq_folder,
                         ref_folder = ref_folder,
                         include_introns = include_introns,
                         slurm_out = paste0(slurm_base, "_count-%j.out"))
    })

    ############ Generate Seurat objects for the things that need Seurat objects
    # Generate Seurat objects and save to SeuratObj folder
    #!!!!!!!!!!!!!!! Should I save an object for each species or just one object?
    to_make_sobj <-
        dplyr::filter(sample_data,
                      make_sobj == TRUE)

    parallel::mclapply(unique(to_make_sobj$Sample_ID),
                       mc.cores = proc_threads,
                       function(s_id) {
        message("Generating Seurat object for ", to_make_sobj$Sample_ID[i], ".")


    })

    # for (sample_name in sub_sample_data$Sample_ID) {
    #     s_obj <-
    #         tenx_load_qc(paste0(counts_folder, "/",
    #                             sample_name,
    #                             "/filtered_feature_bc_matrix"),
    #                      violin_plot = FALSE) %>%
    #         auto_subset()
    # }

    #################################
    ### Need to add code to chmod files

    #################################
    ### Clean up un-needed files

    message("Done!")
}

#' Retrieve sequencing data and store it in
#'
#' @param link_folder Folder containing raw data
#' @param dest_folder Location to store files
#' @param domain Where the raw files are hosted
#' @param user Username to pass to smbclient
#' @param user_group Group the user belongs to
#' @param .pw internal use only
#'
#' @export
#'
#' @examples
#' \dontrun{
#' get_raw_data(link_folder = "210527_Roberts_GSL-AG-2157",
#'              dest_folder = "~/tempDir")
#' }
get_raw_data <- function(link_folder,
                         dest_folder,
                         domain = "//igmdata/igm_roberts",
                         user = Sys.info()[["user"]],
                         user_group = "research",
                         .pw) {

  if (!dir.exists(dest_folder)) {
    system(paste("mkdir", dest_folder))
  }

  # Get list of tar files to use for untaring and to return
  # Eventually I need to avoid entering the password twice, but that is a
  # problem for Future Matt. - good luck, Past Matt
  tar_list <- system(paste("cd ", dest_folder, " ; ",
                           "export LD_LIBRARY_PATH=\"\"; ",
                           "smbclient ",
                           domain, " ",
                           "-U ", user,
                           "%",
                           .pw,
                           " ",
                           "-W ", user_group, " ",
                           "-c ",
                           "'cd ", link_folder, "; ",
                           "ls *.tar'",
                           sep = ""),
                     intern = TRUE) %>%
    stringr::str_subset(".tar") %>%
    stringr::str_replace(" +", "") %>%
    stringr::str_remove(" .+")

  # Get raw data from IGM
  return_val <- system(paste("cd ", dest_folder, " ; ",
                             "export LD_LIBRARY_PATH=\"\"; ",
                             "smbclient ",
                             domain, " ",
                              "-U ", user,
                             "%",
                             .pw,
                             " ",
                             "-W ", user_group, " ",
                             "-c ",
                             "'cd ", link_folder, "; ",
                             "mask \"\"; ",
                             "recurse OFF; ",
                             "prompt OFF; ",
                             "mget *.tar *.md5'",
                             sep = ""))

  if (return_val != 0) {
    stop("Data retrieval failed. Error code ", return_val)
  }

  # Check md5 sums to see if data copied properly.
  check_tar_md5(dest_folder)

  # Untar the downloaded data for further use
  untar_cmd <- paste0("cd ", dest_folder, " ; ",
                      "tar ",
                      "--checkpoint=100000 ",
                      "--checkpoint-action=\"echo= %T %t\" ",
                      "-xf ",
                      stringr::str_c(dest_folder,
                                     tar_list,
                                     sep = "/",
                                     collapse = " "))

  return_val <- system(untar_cmd)

  if (return_val != 0) {
    stop("Untar failed. Error code ", return_val)
  }
  return(tar_list)
}

#' Check that downloaded files match the expected md5sums
#'
#' @param folder Folder containing .tar and .tar.md5 files
#'
#' @keywords internal
#'
check_tar_md5 <- function(folder) {
  message("Checking file md5sums.")
# Make this use srun
  md5_cmd <- paste("md5sum ",
                   folder, "/*.tar",
                   sep = "")

  # Make this into a named vector for ease of comparison
  calc_md5 <- system(md5_cmd, intern = TRUE) %>%
    tibble::as_tibble() %>%
    tidyr::separate(col = "value", sep = "[ ]+", into = c("md5", "file")) %>%
    dplyr::mutate(file = stringr::str_remove(file, ".+/")) %>%
    dplyr::select(file, md5) %>%
    tibble::deframe()

  md5_files <- list.files(path = folder, pattern = ".md5", full.names = TRUE)

  for (md5 in md5_files) {
    md5_true <- readr::read_table2(md5,
                                   col_names = c("md5", "file"),
                                   col_types = "cc")

    if (md5_true$md5 != calc_md5[[md5_true$file]]) {
      stop("md5 checksum check failed for file md5_true$file.")
    }
  }
  message("md5 checksums good!")
}

#' Fix sampleSheet.csv for cellranger mkfastq
#'
#' @param sample_sheet SampleSheet.csv from the BCL folder
#' @param sample_info_sheet Tab delimited sheet with columns containing info to
#'     insert into the sample sheet
#'
#' @keywords internal
#'
fix_sample_sheet <- function(sample_sheet, sample_info_sheet) {
  package_dir <- find.package("rrrSingleCellUtils")

  temp_file <- tempfile(fileext = ".csv")

  system_cmd <- paste0("perl ",
                      package_dir,
                      "/exec/fixSampleSheet.pl ",
                      "--sampleInfo ", sample_info_sheet, " ",
                      "--sampleSheet ", sample_sheet, " ",
                      "> ", temp_file, "; ",
                      "mv ", temp_file, " ", sample_sheet)

  return_val <- return(system(system_cmd))

  if (return_val != 0) {
    stop("Sample sheet repair failed. Error code ", return_val)
  }
}

#' Run cellranger mkfastq
#'
#' @param sample_info File containing sample info (see details)
#' @param tar_folders Character vector containing folders extracted from tar
#' @param bcl_folder Path to write BCL files
#' @param fastq_folder Path to write fastq files
#' @param email Email for Slurm notifications
#' @param slurm_out Location to write out slurm out files
#'
#' @export
#'
#' @examples
#' \dontrun{
#' need example
#' }
cellranger_mkfastq <- function(sample_info,
                               email = "",
                               tar_folders,
                               bcl_folder = "/home/gdrobertslab/lab/BCLs",
                               fastq_folder = "/home/gdrobertslab/lab/FASTQs",
                               slurm_out = paste(getwd(),
                                                 "/slurmOut_mkfastq-%j.out",
                                                 sep = "")) {
  sample_data <- readr::read_delim(sample_info,
                                   delim = "\t",
                                   col_names = TRUE,
                                   trim_ws = TRUE,
                                   show_col_types = FALSE)

  run_name <- sample_data$Sample_Project[1]

  if (email != "") {
    email <- paste("#SBATCH --mail-user=", email, "\n",
                   "#SBATCH --mail-type=ALL\n",
                   sep = "")
  }

  # This assumes that the sampleInfoSheet has a single exp_type.
  # process_raw_data writes out temp sampleInfoSheets broken up by exp_type
  if (sample_data$exp_type %>% unique %>% length() > 1) {
      warning("cellranger_mkfastq() requires a single exp_type per run.")
      warning("Break up sampleInfoSheet by exp_type and try again.")
      stop()
  }

  cellranger_type <- "cellranger"
  if (grepl("multiomics", sample_data$exp_type[1])) {
      cellranger_type <- "cellranger-arc"
  }

  filter_arg <- ""
  if (system(paste0("grep 'PlannedIndex2ReadCycles>0' ",
             bcl_folder, "/",
             sample_data$Sample_Project[1], "/",
             tar_folders[1],
             "/RunParameters.xml > /dev/null ")) == 0) {
    filter_arg <- "\\\\\\\n  --force-single-index\\\n"
  }

  # Need to add a suffix to the fastq folder for the multiomics so we can
  # reference it during counting
    fastq_suffix <- ""
    base_mask <- "\n  --use-bases-mask=Y28n*,I10n*,I10n*,Y90n* \\\\\\\n"
    if (sample_data$exp_type[1] == "multiomics GEX") {
        fastq_suffix <- "_R"
        base_mask <-
        "\n  --use-bases-mask=Y28n*,I10n*,I10n*,Y90n* \\\\\\\n  --filter-dual-index \\\\"
    } else if (sample_data$exp_type[1] == "multiomics ATAC") {
        fastq_suffix <- "_A"
        base_mask <-
            "\n  --use-bases-mask=Y50n*,I8n*,Y24n*,Y49n* \\\\\\\n  --filter-single-index \\\\"
    }

  replace_tibble <-
    tibble::tibble(find = c("placeholder_run_name",
                            "placeholder_array_max",
                            "placeholder_bcl_folder_array",
                            "placeholder_bcl_path",
                            "placeholder_fastq_folder",
                            "placeholder_email",
                            "placeholder_slurm_out",
                            "placeholder_slurm_name",
                            "placeholder_cellranger_type",
                            "placeholder_exp_type",
                            "placeholder_filter_arg",
                            "placeholder_base_mask"),
                   replace = c(run_name,
                               length(tar_folders) - 1,
                               paste(tar_folders, collapse = " "),
                               bcl_folder,
                               fastq_folder,
                               email,
                               slurm_out,
                               paste("mkfastq_",
                                     run_name,
                                     sep = ""),
                               cellranger_type,
                               fastq_suffix,
                               filter_arg,
                               base_mask))

  package_dir <- find.package("rrrSingleCellUtils")

  sbatch_template <-
    readr::read_file(paste(package_dir,
                           "/cellranger_demux_template.job",
                           sep = ""))

  # Replace placeholders with real data
  for (i in seq_len(nrow(replace_tibble))) {
    sbatch_template <-
      stringr::str_replace_all(sbatch_template,
                               pattern = replace_tibble$find[i],
                               replacement = replace_tibble$replace[i])
  }

  temp_file <- tempfile(fileext = ".sh",
                        pattern = "mkfastq",
                        tmpdir = getwd())

  readr::write_file(sbatch_template, file = temp_file)

  return_val <- system(paste("sbatch", temp_file))

  if (return_val != 0) {
    stop("Cellranger mkfastq sbatch submission failed. Error code ", return_val)
  }
}

#' Run cellranger count on 10X data
#'
#' @param sample_info File with sample information Required columns:
#'     Sample_Project, Sample_ID, Reference, Cell_Num
#' @param email Email for Slurm notifications
#' @param include_introns Should intronic reads be included in counts?
#' @param realign_suffix Suffix to add to output folders when re-aligning
#' @param counts_folder Folder for cellranger counts output
#' @param fastq_folder Path to write fastq files
#' @param ref_folder Path to 10x reference folders
#' @param slurm_out Location to write out slurm out files
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Placeholder
#' }
cellranger_count <- function(sample_info,
                             email = "",
                             include_introns = FALSE,
                             realign_suffix = "",
                             counts_folder = "/home/gdrobertslab/lab/Counts",
                             fastq_folder = "/home/gdrobertslab/lab/FASTQs",
                             ref_folder = "/home/gdrobertslab/lab/GenRef",
                             slurm_out = paste(getwd(),
                                               "/slurmOut_count-%j.out",
                                               sep = "")) {

  sample_data <- readr::read_delim(sample_info,
                                   delim = "\t",
                                   col_names = TRUE,
                                   trim_ws = TRUE,
                                   show_col_types = FALSE)

  run_name <- sample_data$Sample_Project[1]

  fastq_folder <- paste(fastq_folder,
                        run_name,
                        sep = "/")

  if (email != "") {
    email <- paste0("#SBATCH --mail-user=", email, "\n",
                    "#SBATCH --mail-type=ALL\n")
  }

  multiomic <- dplyr::if_else(grepl("multiomic", sample_data$exp_type[1]),
                              TRUE,
                              FALSE)

  cellranger_template <- "/cellranger_count_template.job"
  tmp_csv <- ""
  # Need to deal with the multiomic completely differently due to cellranger-arc
  if (multiomic) {
    cellranger_template <- "/cellranger_count_multiomics_template.job"

    tmp_csv <- tempfile(pattern = "crCount",
                        tmpdir = getwd(),
                        fileext = ".csv")

    # Create csv that cellranger-arc can use for count - each sample will be
    # pulled out of this csv in the job template using grep
    sample_data %>%
        dplyr::select(Sample_ID, exp_type) %>%
        dplyr::mutate(suffix = dplyr::if_else(exp_type == "multiomics GEX",
                                              "_R",
                                              "_A"),
                      fastqs = paste0(fastq_folder,
                                      suffix,
                                      "/",
                                      run_name),
                      exp_type = stringr::str_replace(exp_type,
                                                      "multiomics GEX",
                                                      "Gene Expression") %>%
                            stringr::str_replace("multiomics ATAC",
                                                 "Chromatin Accessibility")) %>%
        dplyr::rename(library_type = exp_type,
                      sample = Sample_ID) %>%
        dplyr::select(fastqs, sample, library_type, -suffix) %>%
        readr::write_csv(file = tmp_csv)
  }

  # Since cellranger-arc count pulls in both GEX and ATAC at the same time
  # we don't need to run both, so only need the data from one of the two exp
  # to pull in both, that's why I use only gex data in the replace_tibble below
  if (multiomic) {
    sample_data <-
        sample_data %>%
        dplyr::filter(exp_type == "multiomics GEX")
  }

  # the command for including introns is different for cellranger-arc
  intron_arg <-
    dplyr::if_else(include_introns,
                   "",
                   dplyr::if_else(multiomic,
                                  "\n --gex-exclude-introns \\\\\\",
                                  "\n --include-introns false \\\\\\"))

  # Need very different arguments for cellranger/cellranger-arc
  # Using separate template
  # The Cell Ranger ARC pipeline can only analyze Gene Expression and ATAC data
  # together and the input is a csv file
  replace_tibble <- tibble::tibble(find = c("placeholder_run_name",
                                            "placeholder_array_max",
                                            "placeholder_sample_array_list",
                                            "placeholder_outdir_array_list",
                                            "placeholder_reference_array_list",
                                            "placeholder_num_cells_list",
                                            "placeholder_email",
                                            "placeholder_slurm_out",
                                            "placeholder_reference_folder",
                                            "placeholder_counts_folder",
                                            "placeholder_fastq_folder",
                                            "placeholder_slurm_name",
                                            "placeholder_library_csv",
                                            "placeholder_include_introns"),
                                   replace = c(run_name,
                                               nrow(sample_data) - 1,
                                               paste(sample_data$Sample_ID,
                                                     collapse = " "),
                                               paste(paste0(sample_data$Sample_ID,
                                                            realign_suffix),
                                                     collapse = " "),
                                               paste(sample_data$Reference,
                                                     collapse = " "),
                                               paste(sample_data$Cell_Num,
                                                     collapse = " "),
                                               email,
                                               slurm_out,
                                               ref_folder,
                                               counts_folder,
                                               fastq_folder,
                                               paste("count_",
                                                     run_name,
                                                     sep = ""),
                                               tmp_csv,
                                               intron_arg))

  package_dir <- find.package("rrrSingleCellUtils")

  sbatch_template <-
    readr::read_file(paste0(package_dir,
                            cellranger_template))

  # Replace placeholders with real data
  for (i in seq_len(nrow(replace_tibble))) {
    sbatch_template <-
      stringr::str_replace_all(sbatch_template,
                               pattern = replace_tibble$find[i],
                               replacement = replace_tibble$replace[i])
  }

  temp_file <- tempfile(fileext = ".sh",
                        pattern = "count",
                        tmpdir = getwd())

  readr::write_file(sbatch_template, file = temp_file)

  return_val <- system(paste("sbatch", temp_file))

  if (return_val != 0) {
    stop("Cellranger count sbatch submission failed. Error code ", return_val)
  }
}

#' Add columns to sample_info to show if data are downloaded, processed, etc.
#'
#' @param sample_info File with sample information
#' @param bcl_folder Path for BCL files
#' @param fastq_folder Path for fastq files
#' @param counts_folder Path for counts files
#' @param sobj_folder Path for Seurat objects
#'
#' @keywords internal
#'
#' @return sample_info with additional columns
add_data_status <- function(sample_info,
                            bcl_folder = "/home/gdrobertslab/lab/BCLs",
                            fastq_folder = "/home/gdrobertslab/lab/FASTQs",
                            counts_folder = "/home/gdrobertslab/lab/Counts",
                            sobj_folder = "/home/gdrobertslab/lab/SeuratObj") {
    # Check if Seurat object exists
    sample_info$make_sobj <-
        !file.exists(paste0(sobj_folder, "/",
                            sample_info$Sample_ID,
                            ".qs"))

    # Check if counts data exists
    # Need this to be TRUE even if make_sobj is FALSE
    sample_info$run_cellranger_count <-
        !file.exists(paste0(counts_folder, "/",
                           sample_info$Sample_ID, "/",
                           "filtered_feature_bc_matrix/barcodes.tsv.gz"))

    fastq_folder_suffix <-
        stringr::str_replace_all(sample_info$Protocol,
                                 c("3GEX"    = "",
                                   "CNV"     = "",
                                   "^ATAC$"  = "",
                                   "CITE"    = "",
                                   "^MATAC$" = "_A",
                                   "MGEX"    = "_R"))


    # Check if fastq data exists or if run_cellranger_count is TRUE
    sample_info$run_cellranger_mkfastq <-
        lapply(paste0(fastq_folder, "/",
                      sample_info$Sample_Project,
                      fastq_folder_suffix, "/",
                      sample_info$Sample_Project, "/",
                      sample_info$Sample_ID, "/",
                      sample_info$Sample_ID,
                      "*R1*fastq.gz"),
               function(x) length(Sys.glob(x)) == 0) %>%
        unlist() &
        sample_info$run_cellranger_count

    # Check if bcl data exists or if cellranger_mkfastq is TRUE
    sample_info$download_data <-
        (dir.exists(paste0(bcl_folder, "/",
                          sample_info$Sample_Project, "/",
                          sample_info$link_folder)) |
        sample_info$run_cellranger_mkfastq) &
        !is.na(sample_info$link_folder)

    return(sample_info)
}

#' Get the experiment type from a vector of exp_types
#'
#' Input can be one of "MATAC" "MGEX"  "3GEX" or "MATAC" and "MGEX"
#' Need to return "ATAC" "GEX" "GEX" or "ATAC+GEX" if both "MATAC" and "MGEX"
#'
#' @param x Vector of exp_types
#'
#' @keywords internal
get_exp_type <- function(x) {
    if (all(x == "3GEX")) {
        return("GEX")
    } else if (all(sort(x) == c("MATAC", "MGEX"))) {
        return("ATAC+GEX")
    } else if (all(x == "MATAC")) {
        return("ATAC")
    } else if (all(x == "MGEX")) {
        return("GEX")
    } else {
        stop("Unknown exp_type.")
    }
}

#' Make a Seurat object from 10X data
#'
#' @param s_id Sample ID
#' @param sample_info Sample information
#' @param counts_folder Path to counts folder
#' @param sobj_folder Path to Seurat object folder
#' @param cutoff_hist_folder Path to cutoff histogram folder
#'
#' @keywords internal
make_sobj <- function(s_id = s_id,
                      sample_info = to_make_sobj,
                      counts_folder = counts_folder,
                      sobj_folder = sobj_folder,
                      cutoff_hist_folder = cutoff_hist_folder) {

    sample_data <-
        dplyr::filter(sample_info,
                      Sample_ID == s_id)

    exp_type <- get_exp_type(sample_data$Protocol)

    s_obj <-
        tenx_load_qc(h5_file = paste0(counts_folder, "/",
                                        to_make_sobj$Sample_ID[i],
                                        "/filtered_feature_bc_matrix.h5"),
                        violin_plot = FALSE,
                        exp_type = exp_type)

    # Add metadata from sample_data to s_obj
    # For multiomics where there are two rows, use unique and paste to keep one
    # copy of each value
    for (colname in colnames(sample_data)) {
        sobj[[colname]] <-
            sample_data[[colname]] %>%
            unique() %>%
            paste(collapse = ", ")
    }

    # If GEX present, just process it like normal
    # The process for 3' GEX and multiomics is the same
    if (grepl("GEX", exp_type)) {
        s_obj <-
            process_sobj_gex(s_obj,
                             sample_data,
                             cutoff_hist_folder)
    }
    # If ATAC data, need to do some extra stuff
    if (grepl("ATAC", exp_type)) {
        # subset data

        # process data
    } else if (exp_type == "ATAC+GEX") {
        # subset data

        # process data
    }
}

#' Take a GEX seurat object, subset it down based on cutoffs and process it
#'
#' @param s_obj Seurat object
#' @param sample_data Sample information
#' @param cutoff_hist_folder Path to cutoff histogram folder
#'
#' @keywords internal
#'
#' @return Seurat object
process_sobj_gex <- function(s_obj,
                             sample_data,
                             cutoff_hist_folder) {
    Seurat::DefaultAssay(s_obj) <- "RNA"
    # numbers to use for subsetting the data down
    # My assumption here is that this tibble will have one row at this point
    # since it's a single sample and just 3GEX or MGEX
    subset_table <-
        sample_data %>%
        dplyr::filter(grepl("GEX", Protocol)) %>%
        dplyr::select(dplyr::starts_with("subset_")) %>%
        dplyr::rename_all(~stringr::str_remove(., "subset_"))

    if (all(is.na(subset_table))) {
        # No cutoffs values provided, so autocalculate them
        png(paste0(cutoff_hist_folder, "/",
                    sample_data$Sample_ID[1],
                    "_cutoff_hist.png"))
        s_obj <- auto_subset(s_obj,
                             features = c("nFeature_RNA",
                                          "nCount_RNA",
                                          "percent.mt"))
        dev.off()
    } else {
        # Kick out columns with all NAs
        subset_table <-
            subset_table %>%
            dplyr::select_if(~!all(is.na(.)))

        cutoff_table <-
            subset_table %>%
            tidyr::pivot_longer(names_to = "feature",
                                values_to = "value",
                                cols = dplyr::everything()) %>%
            dplyr::mutate(direction = stringr::str_remove(feature, ".+_") %>%
                          paste0("_val"),
                          feature = stringr::str_remove(feature, "_m[ai][nx]$")) %>%
            tidyr::pivot_wider(names_from = direction,
                               values_from = value)


        png(paste0(cutoff_hist_folder, "/",
                   sample_data$Sample_ID[1],
                   "_cutoff_hist.png"))
        print(feature_hist(s_obj,
                     features = cutoff_table$feature,
                     cutoff_table = cutoff_table))
        dev.off()

        for (column in colnames(subset_table)) {
            descriptor <- stringr::str_remove(column, "_m[ai][nx]$")
            direction <- stringr::str_remove(column, "^.+_")
            cutoff_value <- subset_table[[column]][1]
            if (direction == "min") {
                # subset can't accept a variable to name the column for
                # subsetting, so we use this approach instead
                s_obj <-
                    s_obj[, which(Seurat::FetchData(s_obj, vars = descriptor) >= cutoff_value)]
            } else if (direction == "max") {
                s_obj <-
                    s_obj[, which(Seurat::FetchData(s_obj, vars = descriptor) <= cutoff_value)]
            } else {
                print(subset_table)
                stop("Unknown direction in subset table for ",
                        descriptor, ".")
            }
        }
    }

    # process data
    s_obj <- process_seurat(s_obj)

    return(s_obj)
}



process_sobj_atac <- function(s_obj,
                              sample_data,
                              cutoff_hist_folder,
                              gtf,
                              nucl_cutoff = 4,
                              tss_cutoff = 2,
                              frag_files,
                              default_assay = "ATAC") {
    Seurat::DefaultAssay(s_obj) <- default_assay
    # numbers to use for subsetting the data down
    # My assumption here is that this tibble will have one row at this point
    # since it's a single sample and just MATAC
    subset_table <-
        sample_data %>%
        dplyr::filter(Protocol == "MATAC") %>%
        dplyr::select(dplyr::starts_with("subset_")) %>%
        dplyr::rename_all(~stringr::str_remove(., "subset_"))

    # Add in ATAC specific metadata columns
    s_obj <-
        add_atac_metadata(s_obj,
                          gtf = gtf,
                          nucl_cutoff = nucl_cutoff,
                          tss_cutoff = tss_cutoff,
                          frag_files = frag_files)

    if (all(is.na(subset_table))) {
        sobj <- auto_subset(s_obj,
                            features = c("nFeature_ATAC",
                                         "nCount_ATAC"))

    }
    # process data
    s_obj <- process_seurat_atac(s_obj)

    return(s_obj)
}