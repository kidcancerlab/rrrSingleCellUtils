# Roberts Lab R Single Cell Analysis Utility functions

> These are utility functions for use by the Roberts lab to analyze single-cell and genomics data
> sequencing data.

---

## Functions

### Working and tested functions

  * tenx_load_qc()
    * Load 10X data
    * Formerly tenXLoadQC()
      * Now requires mt_pattern and species_pattern instead of spec
  * gen_cellecta_bc_data()
    * Extract cellecta barcode information from a sam or bam file
  * plot_complex_heatmap()
    * Make a complex heatmap showing nichenet output from FindLigands()
  * find_ligands()
    * Formerly findLigands()
    * Run nichenetr to find ligands potentially inducing receptor-driven gene expression changes
  * find_tar_genes()
    * Formerly findTarGenes()
    * Create a gene list containing putative targets of ligand activity
  * kill_cc()
    * Formerly killCC()
    * Regress out cell cycle effects
  * plot_cc()
    * Plot the proportion of cells in each cell cycle stage

### Still in testing

  * process_ltbc()

### To do

  * Make arguments inherit documentation from other functions instead of re-typing them
  * Add argument to provide path to cellranger
  * Vignette
  * For multiomics, make sure both atac and gex are present in sample input sheet, and not as "sample1_atac" and "sample1_gex", which will fail. Should be "sample1" and "sample1" for both.
  * Try out snp calling: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1863-4
  * In tenx_load_qc, check if all data is filtered out or sobj is empty
  * Figure out multiomic FRiP calculation and see if I need to divide by 2
    * https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/algorithms/overview
  * Make optimize_silhouette warn if there are no dim reductions or clustering present
  * Try out snp calling: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1863-4
  * In tenx_load_qc, check if all data is filtered out or sobj is empty

---

## Errors
> If you run into any errors, please run the following commands and send the
> output to me along with the error messages produced.
>
> ```
> traceback()
> sessionInfo()
> ```
>
