#' Download and pre-process raw 10X data
#'
#' @param sample_info File containing sample info (see details)
#' @param link_folder Folder containing raw data
#' @param domain Where the raw files are hosted
#' @param exp_type Type of 10X sequence data
#' @param bcl_folder Path to write BCL files
#' @param fastq_folder Path to write fastq files
#' @param counts_folder Path to write counts files
#'
#' @details Need to put this in
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' placeholder
#' }
process_raw_data <- function(sample_info,
                             link_folder,
                             domain = "//igmdata/igm_roberts",
                             exp_type = "scRNA-seq_3prime",
                             bcl_folder = "/home/gdrobertslab/lab/BCLs",
                             fastq_folder = "/home/gdrobertslab/lab/FASTQs",
                             counts_folder = "/home/gdrobertslab/lab/Counts") {

  #project, from sample_info file?
  dest_folder <- paste(bcl_folder, "/",
                       project,
                       sep = "")

  # Download, untar and check md5 sums of raw data
  get_raw_data(link_folder = link_folder, dest_folder = dest_folder)

  # Fix sample sheet for all folders with BCL files
  list.files(path = dest_folder, pattern = ".tar$") %>%
    stringr::str_remove(".tar") %>%
    lapply(., function(x)
      fix_sample_sheet(sample_sheet = paste(dest_folder, "/",
                                            x, "/SampleSheet.csv",
                                            sep = ""),
                       sample_info_sheet = sample_info))

  # Run cellranger mkfastq


  # Run cellranger count if scRNA-seq_3prime
  cellranger_count(sample_info = sample_info)

  # Run cellranger-DNA if
}

#' Retrieve sequencing data and store it in
#'
#' @param link_folder Folder containing raw data
#' @param dest_folder Location to store files
#' @param domain Where the raw files are hosted
#' @param user Username to pass to smbclient
#' @param user_group Group the user belongs to
#'
#' @return
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

  get_data_cmd <- paste("cd ", dest_folder, " ; ",
                        "export LD_LIBRARY_PATH=\"\"; ",
                        "smbclient ",
                        domain, " ",
                        "-U ", user,
                        "%", getPass::getPass("Password for smbclient"), " ",
                        "-W ", user_group, " ",
                        "-c ",
                        "'cd ", link_folder, "; ",
                        "mask \"\"; ",
                        "recurse OFF; ",
                        "prompt OFF; ",
                        "mget *.tar *.md5'",
                        sep = "")

  system(get_data_cmd)

  check_tar_md5(dest_folder)

  untar_cmd <- paste("cd ", dest_folder, " ; ",
                     "tar ",
                     "--checkpoint=100000 ",
                     "--checkpoint-action=\"echo= %T %t\" ",
                     "-xf ",
                     dest_folder,
                     "/*tar",
                     sep = "")

  system(untar_cmd)
}

#' Check that downloaded files match the expected md5sums
#'
#' @param folder Folder containing .tar and .tar.md5 files
#'
#' @keywords internal
#'
check_tar_md5 <- function(folder) {
  message("Checking file md5sums.")

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

  system_cmd <- paste(package_dir,
                      "/exec/fixSampleSheet.pl ",
                      "--sampleInfo ", sample_info_sheet,
                      "--sampleSheet ", sample_sheet,
                      " > ", sample_sheet,
                      sep = "")

  system(system_cmd)
}

#' Run cellranger count on 10X data
#'
#' @param sample_info File with sample information Required columns:
#'     Sample_Project, Sample_ID, Reference, Cell_Num
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Placeholder
#' }
cellranger_count <- function(sample_info) {
  sbatch_template <- readr::read_file("inst/cellranger_count_template.job")

  sample_data <- readr::read_delim(sample_info,
                                   delim = "\t",
                                   col_names = TRUE,
                                   trim_ws = TRUE)

  replace_tibble <- tibble::tibble(find = c("placeholder_run_name",
                                    "placeholder_array_max",
                                    "placeholder_sample_array_list",
                                    "placeholder_reference_array_list",
                                    "placeholder_num_cells_list"),
                           replace = c(unique(sample_data$Sample_Project),
                                       nrow(sample_data) %>%
                                         as.character(),
                                       paste(sample_data$Sample_ID,
                                             collapse = " "),
                                       paste(sample_data$Reference,
                                             collapse = " "),
                                       paste(sample_data$Cell_Num,
                                             collapse = " "))
  )

  # Replace placeholders with real data
  for (i in seq_len(nrow(replace_tibble))) {
    sbatch_template <-
      stringr::str_replace_all(sbatch_template,
                               pattern = replace_tibble$find[i],
                               replacement = replace_tibble$replace[i])
  }

  temp_file <- tempfile(fileext = ".sh")

  readr::write_file(sbatch_template, file = temp_file)

  system(paste("sbatch", temp_file))
}
