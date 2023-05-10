# Roberts Lab R Single Cell Analysis Utility functions

> These are utility functions for use by the Roberts lab to analyze single cell
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
  * process_raw_data
    * This is a wrapper for several functions to get and process single
    cell data from the NCH IGM core. I have built the defaults to be specific to
    the Roberts lab, so you may need to carefully change the defaults if you
    want to use it outside of this context. The input is a link to data and a
    sample sheet that outlines the information about each sample (see
    inst/exampleSampleInfoSheet.txt - Column headers must remain unchanged). The
    data are then downloaded using smbclient and then md5sum checked and
    untar’d. The data are then processed with cellranger mkfastq and either
    cellranger count or cellranger-dna cnv (depending on the exp_type argument).

### Still in testing

  * process_ltbc()

### To do

  * Make arguments inherit from other functions instead of re-typing them
  * Add argument to provide path to cellranger
  * Add argument to specify bcl2fastq path or module
  * Vignette
  * Function to check sample sheet
  * For multiomics, make sure both atac and gex are present in sample input sheet, and not as "sample1_atac" and "sample1_gex", which will fail. Should be "sample1" and "sample2" for both.
  * Save shell scripts to current directory instead of tempdir

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
