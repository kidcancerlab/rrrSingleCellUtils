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
  
### Still in testing

  * process_ltbc()

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
