#' Create a gene list containing putative targets of ligand activity.
#'
#' @details This gene list is needed as one of the inputs for the
#'    ligand-receptor analysis. This list should be generated to identify genes
#'    in the same target cells used in the findLigands function as the "targets"
#'    parameter. You do not have to use this particular function to generate
#'    that gene list, but you can use this function to do it.
#' @param sobject A Seurat object containing all of the cells for analysis
#'    (required)
#' @param id1 The idents of the target cells (ligand-stimulated cells)
#'    (required)
#' @param id2 The idents of the unstimlated cells for comparison (required)
#' @param pval The p-value to use as a cutoff for up-regulation of genes
#'    (default = 0.05)
#' @param logfc The log2 fold-change to use as cutoff for up-regulation of
#'    genes (default = 0.25)
#' @param spec The species of the gene set (default = "human", can also be
#'    "mouse")
#' @export
#'
#' @return A list of target genes
#'
#' @examples
#' \dontrun{
#' targets <- find_tar_genes(seurat_obj,
#'                           id1 = "d14",
#'                           id2 = "d35",
#'                           logfc = 0.25,
#'                           spec = "human")
#' }

find_tar_genes <- function(sobject, id1, id2, pval = 0.05, logfc = 0.25,
                         spec = "human") {

  if (spec == "human") {
    ligand_list <- human_ligand_list
  } else if (spec == "mouse") {
    ligand_list <- mouse_ligand_list
  } else {
    stop("for findTarGenes, spec must be defined as either human or mouse")
  }

  targets <-
    Seurat::FindMarkers(sobject, ident.1 = id1, ident.2 = id2) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(p_val_adj <= pval & abs(avg_log2FC) >= logfc) %>%
    dplyr::pull(gene) %>%
    dplyr::intersect(ligand_list)

  return(targets)
}


#' Function findLigands
#'
#' @param sobject
#' @param gset
#' @param receiver
#' @param senders
#' @param gset_spec
#' @param rec_pct
#' @param rec_spec
#' @param send_pct
#' @param send_spec
#' @param stringency
#' @param dPlot
#' @param LTVis
#' @param LRVis
#' @param nBest
#'
#' @details
#' Identify ligands expressed by sender cells that activate target genes in
#'  receiver cells.
#  This workflow leverages the nicheNET code written by Robin Browaeys and
#' published in https://doi.org/10.1038/s41592-019-0667-5
#' Note that this function will translate all murine genesets into their human
#'  orthologs for processing. Results of analyses using murine samples should
#'  be vetted and evaluated for legitimacy.
#'
#' @return A list containing the following objects: dotplots for
#'   1. ligand dotplot
#'   2. receptor dotplot (ggplot, if selected)
#'   3. ligand-target graph (ggplot, if selected)
#'   4. ligand-receptor interaction graph (ggplot, if selected)
#'   5. ligand activity matrix (tibble)
#'   6. ligand-target matrix for the prioritized ligands (tibble)
#'   7. ligand-receptor interaction matrix (tibble)
#'
#' @export
#'
#' @examples
#' \donrun{
#' output <- find_ligands(combo, gset = targets_2,receiver = "Anchors",
#'                        senders = c("M1-like",
#'                                    "M2-like",
#'                                    "AlvMac",
#'                                    "DC"),
#'                        gset_spec = "human",
#'                        rec_spec = "human",
#'                        send_spec = "mouse")
#' }
find_ligands <- function(sobject, gset, receiver, senders, gset_spec = "human",
                        rec_pct = 0.10, rec_spec = "human", send_pct = 0.10,
                        send_spec = "human", stringency = "loose", dPlot = TRUE,
                        LTVis = TRUE, LRVis = TRUE, nBest = 20) {

  # Set placeholders for optional objects
  dp = NULL
  ltv = NULL
  vis_ligand_target = NULL

  # Make a list of expressed genes in the receiver cells, translating mouse to
  # human if needed
  expressed_genes_receiver <- nichenetr::get_expressed_genes(receiver,
                                                             sobject,
                                                             pct = rec_pct)

  if (rec_spec == "mouse") {

    ################################ Make function and putt all calls together
    expressed_genes_receiver <-
      nichenetr::convert_mouse_to_human_symbols(expressed_genes_receiver) %>%
      na.omit () %>%
      as.character()
  }

  background_expressed_genes <- intersect(expressed_genes_receiver,
                                          rownames(ligand_target_matrix))

  # Make a list of expressed genes in the sender cells
  # Can take a single value or a vector of sender cluster idents
  expressed_genes_sender <- senders %>%
    unique() %>%
    lapply(nichenetr::get_expressed_genes, sobject, pct = send_pct) %>%
    unlist() %>%
    unique()

  if (send_spec == "mouse") {
    expressed_genes_sender <-
      nichenetr::convert_mouse_to_human_symbols(expressed_genes_sender) %>%
      na.omit() %>%
      as.character()
  }

  # Identify the ligands and receptors within senders and receivers and find
  # L-R pairs
  # Analysis is sensitive to either loose or strict stringency with regard to
  # the evidence behind the L-R interactions

  ########################## Simplify copied/pasted code
  if (stringency == "loose") {
    expressed_receptors <- intersect(receptors, expressed_genes_receiver)
    expressed_ligands <- intersect(ligands, expressed_genes_sender)
    potential_ligands <- lr_network %>%
      filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
      pull(from) %>%
      unique()
  } else if (stringency == "strict") {
    expressed_receptors <- intersect(receptors_bona_fide, expressed_genes_receiver)
    expressed_ligands <- intersect(ligands_bona_fide, expressed_genes_sender)
    potential_ligands <- lr_network_strict %>%
      filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
      pull(from) %>%
      unique()
  } else {
    stop("This function requires that the parameter `stringency` be set to
         either `loose` or `strict`, the default being `loose`.")
  }

  if (is_empty(potential_ligands)) {
    stop("No potential receptor-ligand pairs are identified in these cells with
         these settings.")
  }

  # Convert target gene set into human, if necessary
  if (gset_spec == "mouse") {
    gset <- convert_mouse_to_human_symbols(gset) %>%
      na.omit() %>%
      as.character()
  }

  # Evaluate ligands based on their potential to activate genes within the
  # target gene set
  ligand_activities <-
    nichenetr::predict_ligand_activities(
      geneset = gset,
      background_expressed_genes = background_expressed_genes,
      ligand_target_matrix = ligand_target_matrix,
      potential_ligands = potential_ligands)

  ligand_activities <- ligand_activities %>%
    dplyr::arrange(-pearson) %>%
    mutate(rank = rank(desc(pearson)))

  best_upstream_ligands <- ligand_activities %>%
    dplyr::top_n(nBest, pearson) %>%
    dplyr::arrange(-pearson) %>%
    dplyr::pull(test_ligand) %>%
    unique()

  # Generate a matrix of the ligands and their downstream gene targets
  active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(nichenetr::get_weighted_ligand_target_links,
           geneset = gset,
           ligand_target_matrix = ligand_target_matrix,
           n = 200) %>%
    bind_rows()

  active_ligand_target_links <-
    nichenetr::prepare_ligand_target_visualization(
      ligand_target_df = active_ligand_target_links_df,
      ligand_target_matrix = ligand_target_matrix,
      cutoff = 0.33)

  order_ligands <- intersect(best_upstream_ligands,
                            colnames(active_ligand_target_links)) %>%
    rev() %>%
    make.names()

  order_targets <- active_ligand_target_links_df$target %>%
    unique() %>%
    intersect(rownames(active_ligand_target_links)) %>%
    make.names()

  rownames(active_ligand_target_links) <-
    rownames(active_ligand_target_links) %>%
    make.names() # make.names() for heatmap visualization of genes like H2-T23

  colnames(active_ligand_target_links) <-
    colnames(active_ligand_target_links) %>%
    make.names() # make.names() for heatmap visualization of genes like H2-T23

  vis_ligand_target = active_ligand_target_links[order_targets,
                                                 order_ligands] %>%
    t()

  if (gset_spec == "mouse") {
    colnames(vis_ligand_target) <-
      convert_human_to_mouse_symbols(colnames(vis_ligand_target))
  }
  if (send_spec == "mouse") {
    rownames(vis_ligand_target) <-
      convert_human_to_mouse_symbols(rownames(vis_ligand_target))
  }

  vis_ligand_target <- vis_ligand_target[!is.na(rownames(vis_ligand_target)),
                                         !is.na(colnames(vis_ligand_target))]

  # If requested, create a ligand-target visualization graph
  if (LTVis == TRUE) {
    ltv <- vis_ligand_target %>%
      make_heatmap_ggplot("Prioritized ligands",
                          "Predicted target genes",
                          color = "blue",
                          legend_position = "top",
                          x_axis_position = "top",
                          legend_title = "Regulatory potential") +
      theme(axis.text.x = element_text(face = "italic")) +
      scale_fill_gradient2(low = "steelblue",
                           high = "steelblue4",
                           breaks = c(0,0.006,0.012))
    print (ltv)
  }

  # Generate a matrix of the predicted ligand-receptor interactions
  lr_network_top = lr_network %>%
    filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
    distinct(from,to)

  best_upstream_receptors = lr_network_top %>%
    pull(to) %>%
    unique()

  lr_network_top_df_large = weighted_networks_lr %>%
    filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

  lr_network_top_df = lr_network_top_df_large %>%
    spread("from","weight",fill = 0)

  lr_network_top_matrix = lr_network_top_df %>%
    select(-to) %>%
    as.matrix() %>%
    magrittr::set_rownames(lr_network_top_df$to)

  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]

  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

  order_receptors = order_receptors %>%
    intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>%
    intersect(colnames(lr_network_top_matrix))

  vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors,
                                                      order_ligands_receptor]
  if (send_spec == "mouse") {
    colnames(vis_ligand_receptor_network) <-
      as.character(convert_human_to_mouse_symbols(colnames(vis_ligand_receptor_network)))
  }
  if (rec_spec == "mouse") {
    rownames(vis_ligand_receptor_network) <-
      as.character(convert_human_to_mouse_symbols(rownames(vis_ligand_receptor_network)))
  }

  vis_ligand_receptor_network <-
    vis_ligand_receptor_network[!is.na(rownames(vis_ligand_receptor_network)),
                                !is.na(colnames(vis_ligand_receptor_network))]

  # If requested, create dot plots of the receptors and ligands
  if (dPlot == TRUE) {
    dotplot_reciever <- DotPlot(subset(sobject, idents = receiver),
                   features = rev(rownames(vis_ligand_receptor_network)),
                   cols = c("gray95", "sienna3")) &
      RotatedAxis()
    print (dotplot_reciever)
    dotplot_ligand <- DotPlot(subset(sobject, idents = senders),
                   features = rev(colnames(vis_ligand_receptor_network)),
                   cols = c("gray95", "sienna3")) &
      RotatedAxis() & coord_flip()
    print (dotplot_ligand)
  }

  # If requested, create a ligand-receptor visualization graph, weighted by
  # transcriptional changes
  if (LRVis == TRUE) {
    lrv <- vis_ligand_receptor_network %>%
      t() %>%
      make_heatmap_ggplot("Ligands",
                          "Receptors",
                          x_axis_position = "top",
                          legend_title = "Prior interaction potential") +
      scale_fill_gradient (low = "white", high = "darkseagreen4")
    print (lrv)
  }

  return.list <- list(dotplot_ligand,
                      dotplot_reciever,
                      ltv,
                      lrv,
                      ligand_activities,
                      vis_ligand_target,
                      vis_ligand_receptor_network)

}
