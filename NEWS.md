# rrrSingleCellUtils 0.19.0
- Updated process_seurat() to include an argument to specify the number of dimensions to use when running RunPCA().
- Suppressed superfluous messages when running annotate_celltypes() from the celldex package.

# rrrSingleCellUtils 0.18.0
- Fixed a bug in merge_atac() where bed entries were not filtered properly.

# rrrSingleCellUtils 0.17.0

- Added an option to find_ligand() to not show the plots that are generated. Also added names to the returned list of results.

# rrrSingleCellUtils 0.16.0

- Added in the option to specify the species (currently either mouse or human) when running `kill_cc()`. Also added the option not to plot the results of the function.

# rrrSingleCellUtils 0.15.0

- Added in functionality and a vignette written by Matt Gust to do velocity analysis on a Seurat object.

# rrrSingleCellUtils 0.14.0

- Removed all of the functions associated with automated processing of raw data. This code is now in the rrrAutoProcess package. This was done as the code was getting too complex and the code was also very specifically tailored to our software environment, which really limits it's general use. This package will now focus on post processing of data and analysis of data.

# rrrSingleCellUtils 0.13.0

- Added in new function annotate_celltypes() which is a wrapper for SingleR to annotate cell types in a Seurat object.

# rrrSingleCellUtils 0.12.0

- Added in the `run_fdl()` function to run force directed layout analysis on a seurat object.
- Also fixed version number since it apparently got mixed up last release.

# rrrSingleCellUtils 0.10.1

- Fixed an issue where filtering based on species_pattern could remove all the genes or peaks in a dataset and result in a rather uninformative error. This has been updated so an informative error is thrown if this happens.

# rrrSingleCellUtils 0.11.1

- Fixed an issue where filtering based on species_pattern could remove all the genes or peaks in a dataset and not result in an error being thrown. As a result a rather uninformative error was given. This has been updated so an informative error is thrown if this happens.

# rrrSingleCellUtils 0.10.0

- Merge in some branches related to SNV calling

# rrrSingleCellUtils 0.9.0

# rrrSingleCellUtils 0.8.0

- Added the option to load in data using an h5 file, which is needed when working with multiomics data.
- I also added another option for frag_file which is again for multiomics
- Did a bit of under the hood work implementing those.
- I also started putting together vignettes to both show how to use all the functions in the package and as something I can run to test that the code is working properly until I get testthat implemented.
