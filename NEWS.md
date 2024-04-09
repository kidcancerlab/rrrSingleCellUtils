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
