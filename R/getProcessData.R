#' Download and pre-process raw 10X data
#'
#' @param sample_info File containing sample info (see details)
#' @param domain Where the raw files are hosted
#' @param email Email for Slurm notifications
#' @param slurm_base Base name for slurm output files
#' @param bcl_folder Path to write BCL files
#' @param fastq_folder Path to write fastq files
#' @param counts_folder Path to write counts files
#' @param ref_folder Path to 10x reference folders
#'
#' @details This is a wrapper for several functions to get and process single
#'  cell data from the NCH IGM core. I have built the defaults to be specific to
#'  the Roberts lab, so you may need to carefully change the defaults if you
#'  want to use it outside of this context. The input is a link to data and a
#'  sample sheet that outlines the information about each sample. The data
#'  are then downloaded using smbclient and then md5sum checked and untar’d.
#'  The data are then processed with cellranger mkfastq and either cellranger
#'  count or cellranger-dna cnv (depending on the exp_type argument).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' process_raw_data(sample_info = "testSampleInfoSheet_test.txt",
#' link_folder = "210913_Roberts_GSL-AG-2350",
#' email = "matthew.cannon@nationwidechildrens.org")
#' }
process_raw_data <- function(sample_info,
                             domain = "//igmdata/igm_roberts",
                             email = "",
                             slurm_base = paste(getwd(), "/slurmOut", sep = ""),
                             bcl_folder = "/home/gdrobertslab/lab/BCLs",
                             fastq_folder = "/home/gdrobertslab/lab/FASTQs",
                             counts_folder = "/home/gdrobertslab/lab/Counts",
                             ref_folder = "/home/gdrobertslab/lab/GenRef") {

  system(paste0("dos2unix ", sample_info))

  sample_data <- readr::read_delim(sample_info,
                                   delim = "\t",
                                   col_names = TRUE,
                                   trim_ws = TRUE,
                                   show_col_types = FALSE)

  project <- sample_data$Sample_Project[1]

  dest_folder <- paste(bcl_folder, "/",
                       project,
                       sep = "")

  exp_type <- sample_data$exp_type[1]

  # Sometimes we're going to have multiple download folders and multiple
  # experiment types (such as multiomics) that need to be handled separately
  # Therefore, we loop over each link_folder
  # This code is going to assume that each download folder has a single exp_type
  temp_sample_info_sheets <- list()

  link_exp_type <- sample_data %>%
    dplyr::select(exp_type, link_folder) %>%
    dplyr::distinct()

#############################
# Going to need to change up the loop here, you need to re-do your password mid-run

  for (i in seq_len(nrow(link_exp_type))) { #unique(sample_data$link_folder)) {
    dl_link <- link_exp_type$link_folder[i]
    # Download, untar and check md5 sums of raw data
    message("Getting raw data from ", dl_link, ".")
    tar_list <- get_raw_data(link_folder = dl_link,
                             dest_folder = dest_folder)

    # Write out subsampled sampleInfoSheet with data from one link_folder
    temp_sample_info_sheets[[dl_link]] <- tempfile(fileext = ".tsv")
    sample_data %>%
        dplyr::filter(link_folder == dl_link) %>%
        readr::write_tsv(file = temp_sample_info_sheets[[dl_link]])

    # Fix sample sheet for all folders with BCL files
    message("Fixing sample sheets in ", dest_folder, ".")
    tar_list %>%
        stringr::str_remove(".tar") %>%
        lapply(., function(x)
            fix_sample_sheet(sample_sheet = paste(dest_folder, "/",
                                                  x, "/SampleSheet.csv",
                                                  sep = ""),
                             sample_info_sheet = temp_sample_info_sheets[[dl_link]]))

    # Run cellranger mkfastq
    message("Submitting slurm command to create fastq",
            "files using cellranger mkfastq.")
    cellranger_mkfastq(sample_info = temp_sample_info_sheets[[dl_link]],
                       email = email,
                       tar_folders = tar_list %>% stringr::str_remove(".tar"),
                       bcl_folder = bcl_folder,
                       fastq_folder = fastq_folder,
                       slurm_out = paste(slurm_base,
                                         "_mkfastq-%j.out",
                                         sep = ""))
  }

  # Run cellranger count if scRNA-seq_3prime
  if (exp_type == "scRNA-seq_3prime") {
    message("Submtting slurm command to run cellranger count.\n",
            "Slurm messages output to", paste(slurm_base,
                                              "_count-%j.out",
                                              sep = ""))
    cellranger_count(sample_info = sample_info,
                     email = email,
                     counts_folder = counts_folder,
                     fastq_folder = fastq_folder,
                     ref_folder = ref_folder,
                     slurm_out = paste(slurm_base, "_count-%j.out", sep = ""))
  } else if (grepl("multiomics", exp_type)) {
    cellranger_count(sample_info = sample_info,
                     email = email,
                     counts_folder = counts_folder,
                     fastq_folder = fastq_folder,
                     ref_folder = ref_folder,
                     slurm_out = paste(slurm_base, "_count-%j.out", sep = ""))
  } else if (exp_type == "scDNA_CNV") {
    warning("Not yet implimented")
  }
  message("Done!")
}

#' Retrieve sequencing data and store it in
#'
#' @param link_folder Folder containing raw data
#' @param dest_folder Location to store files
#' @param domain Where the raw files are hosted
#' @param user Username to pass to smbclient
#' @param user_group Group the user belongs to
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
                         user_group = "research") {

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
                           getPass::getPass("Password for smbclient: "),
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
                             getPass::getPass("Password for smbclient: "),
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

  system_cmd <- paste("perl ",
                      package_dir,
                      "/exec/fixSampleSheet.pl ",
                      "--sampleInfo ", sample_info_sheet, " ",
                      "--sampleSheet ", sample_sheet, " ",
                      "> ", temp_file, "; ",
                      "mv ", temp_file, " ", sample_sheet,
                      sep = "")

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
                            "placeholder_filter_arg"),
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
                               stringr::str_replace_all(sample_data$exp_type[1],
                                                        " ",
                                                        "_"),
                               filter_arg))

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

  temp_file <- tempfile(fileext = ".sh")

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

  cellranger_template <- "/cellranger_count_template.job"
  tmp_csv <- ""
  if (grepl("multiomic", sample_data$exp_type[1])) {
    cellranger_template <- "/cellranger_count_multiomics_template.job"

    tmp_csv <- tempfile(pattern = "crCount",
                        tmpdir = getwd(),
                        fileext = ".csv")

    sample_data %>%
        dplyr::select(Sample_ID, exp_type) %>%
        dplyr::mutate(fastqs = paste0(fastq_folder,
                                      "_",
                                      stringr::str_replace_all(exp_type,
                                                                  " ",
                                                                  "_"),
                                      "/",
                                      run_name),
               exp_type = stringr::str_replace(exp_type,
                                               "multiomics GEX",
                                               "Gene Expression") %>%
                          stringr::str_replace("multiomics ATAC",
                                               "Chromatin Accessibility")) %>%
        dplyr::rename(library_type = exp_type,
                      sample = Sample_ID) %>%
        dplyr::select(fastqs, sample, library_type) %>%
        readr::write_csv(file = tmp_csv)
  }

  # Since cellranger-arc count pulls in both GEX and ATAC at the same time
  # we don't need to run both, so only need the data from one of the two exp
  # to pull in both, that's why I use only gex data in the replace_tibble below
  if (grepl("multiomics", sample_data$exp_type[1])) {
    sample_data <-
        sample_data %>%
        dplyr::filter(exp_type == "multiomics GEX")
  }

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
                                            "placeholder_library_csv"),
                                   replace = c(run_name,
                                               nrow(sample_data) - 1,
                                               paste(sample_data$Sample_ID,
                                                     collapse = " "),
                                               paste(paste(sample_data$Sample_ID,
                                                           realign_suffix,
                                                           sep = ""),
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
                                               tmp_csv))

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

  temp_file <- tempfile(fileext = ".sh")

  readr::write_file(sbatch_template, file = temp_file)

  return_val <- system(paste("sbatch", temp_file))

  if (return_val != 0) {
    stop("Cellranger count sbatch submission failed. Error code ", return_val)
  }
}
