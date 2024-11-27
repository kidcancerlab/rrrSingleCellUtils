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
#' @param sobject A Seurat object containing all of the cells for analysis
#' @param gset A character vector of gene symbols that represent genes altered
#' in the target condition within the receiver cells, can come from findTarGenes
#' function
#' @param receiver The identity of the (single) cluster being stimulated by some
#' ligand
#' @param senders The identities of the cluster(s) producing the ligands
#' @param gset_spec The species of the target gene set
#' @param rec_pct A fraction between 0 and 1, the threshold for inclusion,
#' receptors must be expressed in this fraction of receiver cells
#' @param rec_spec The species of the receiver cells, "human" or "mouse"
#' @param send_pct A fraction between 0 and 1, the threshold for inclusion,
#' ligands must be expressed in this fraction of sender cells
#' @param send_spec The species of the sender cells, "human" or "mouse"
#' @param stringency Determine whether the ligand-receptor interactions are
#' interpreted strictly (based on bona-fide interactions from the literature) or
#' loosely (including predicted and inferred interaction) (must be either
#' "strict" or "loose")
#' @param d_plot Logical, to trigger output of a dotplot of ligands
#' @param lt_vis Logical, to trigger output of a ligand-target visualization
#' graph
#' @param lr_vis Logical, to trigger output of a ligand-receptor visualization
#' graph, weighted by downstream changes in gene expression
#' @param n_best An integer, the number of top ligands to include on the graphs
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
#' \dontrun{
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
                        send_spec = "human", stringency = "loose",
                        d_plot = TRUE, lt_vis = TRUE, lr_vis = TRUE,
                        n_best = 20) {

  if (is.null(rrr_env$ligands)) {
    load_lig_receptor_data()
  }

  # Make a list of expressed genes in the receiver cells, translating mouse to
  # human if needed
  expressed_genes_receiver <- nichenetr::get_expressed_genes(receiver,
                                                             sobject,
                                                             pct = rec_pct)

  if (rec_spec == "mouse") {
    expressed_genes_receiver <-
      nichenetr::convert_mouse_to_human_symbols(expressed_genes_receiver) %>%
      stats::na.omit() %>%
      as.character()
  }

  background_expressed_genes <-
    intersect(expressed_genes_receiver,
              rownames(rrr_env$ligand_target_matrix))

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
      stats::na.omit() %>%
      as.character()
  }

  # Identify the ligands and receptors within senders and receivers and find
  # L-R pairs
  # Analysis is sensitive to either loose or strict stringency with regard to
  # the evidence behind the L-R interactions

  if (stringency == "loose") {
    expressed_receptors <- intersect(rrr_env$receptors,
                                     expressed_genes_receiver)
    expressed_ligands <- intersect(rrr_env$ligands, expressed_genes_sender)
    potential_ligands <- rrr_env$lr_network %>%
      dplyr::filter(from %in% expressed_ligands &
                      to %in% expressed_receptors) %>%
      dplyr::pull(from) %>%
      unique()
  } else if (stringency == "strict") {
    expressed_receptors <- intersect(rrr_env$receptors_bona_fide,
                                     expressed_genes_receiver)
    expressed_ligands <- intersect(rrr_env$ligands_bona_fide,
                                   expressed_genes_sender)
    potential_ligands <- rrr_env$lr_network_strict %>%
      dplyr::filter(from %in% expressed_ligands &
                      to %in% expressed_receptors) %>%
      dplyr::pull(from) %>%
      unique()
  } else {
    stop("This function requires that the parameter `stringency` be set to
         either `loose` or `strict`, the default being `loose`.")
  }

  if (length(potential_ligands) == 0) {
    stop("No potential receptor-ligand pairs are identified in these cells with
         these settings.")
  }

  # Convert target gene set into human, if necessary
  if (gset_spec == "mouse") {
    gset <- nichenetr::convert_mouse_to_human_symbols(gset) %>%
      stats::na.omit() %>%
      as.character()
  }

  # Evaluate ligands based on their potential to activate genes within the
  # target gene set
  ligand_activities <-
    nichenetr::predict_ligand_activities(
      geneset = gset,
      background_expressed_genes = background_expressed_genes,
      ligand_target_matrix = rrr_env$ligand_target_matrix,
      potential_ligands = potential_ligands)

  ligand_activities <- ligand_activities %>%
    dplyr::arrange(-pearson) %>%
    dplyr::mutate(rank = rank(dplyr::desc(pearson)))

  best_upstream_ligands <- ligand_activities %>%
    dplyr::top_n(n_best, pearson) %>%
    dplyr::arrange(-pearson) %>%
    dplyr::pull(test_ligand) %>%
    unique()

  # Generate a matrix of the ligands and their downstream gene targets
  active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(nichenetr::get_weighted_ligand_target_links,
           geneset = gset,
           ligand_target_matrix = rrr_env$ligand_target_matrix,
           n = 200) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(!is.na(weight))

  if (nrow(active_ligand_target_links_df) == 0) {
    stop("No ligand:receptor matches found")
  }

  active_ligand_target_links <-
    nichenetr::prepare_ligand_target_visualization(
      ligand_target_df = active_ligand_target_links_df,
      ligand_target_matrix = rrr_env$ligand_target_matrix,
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

  vis_ligand_target <- t(active_ligand_target_links[order_targets,
                                                    order_ligands])

  if (gset_spec == "mouse") {
    colnames(vis_ligand_target) <-
      nichenetr::convert_human_to_mouse_symbols(order_targets)
  }
  if (send_spec == "mouse") {
    rownames(vis_ligand_target) <-
      nichenetr::convert_human_to_mouse_symbols(order_ligands)
  }

  # by checking first we can avoid messing up a 1 x 1 matrix
  if (NA %in% c(rownames(vis_ligand_target), colnames(vis_ligand_target))) {
    vis_ligand_target <- vis_ligand_target[!is.na(rownames(vis_ligand_target)),
                                           !is.na(colnames(vis_ligand_target))]
  }

  # If requested, create a ligand-target visualization graph
  if (lt_vis == TRUE) {
    ltv <- vis_ligand_target %>%
      nichenetr::make_heatmap_ggplot("Prioritized ligands",
                                     "Predicted target genes",
                                     color = "blue",
                                     legend_position = "top",
                                     x_axis_position = "top",
                                     legend_title = "Regulatory potential") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"))
    print(ltv)
  } else {
    ltv <- NULL
  }

  # Generate a matrix of the predicted ligand-receptor interactions
  lr_network_top <- rrr_env$lr_network %>%
    dplyr::filter(from %in% best_upstream_ligands &
                    to %in% expressed_receptors) %>%
    dplyr::distinct(from, to)

  best_upstream_receptors <- lr_network_top %>%
    dplyr::pull(to) %>%
    unique()

  lr_network_top_df_large <- rrr_env$weighted_networks_lr %>%
    dplyr::filter(from %in% best_upstream_ligands &
                    to %in% best_upstream_receptors)

  lr_network_top_df <- lr_network_top_df_large %>%
    tidyr::spread("from", "weight", fill = 0)

  lr_network_top_matrix <- lr_network_top_df %>%
    dplyr::select(-to) %>%
    as.matrix() %>%
    magrittr::set_rownames(lr_network_top_df$to)

  dist_receptors <- stats::dist(lr_network_top_matrix, method = "binary")
  hclust_receptors <- stats::hclust(dist_receptors, method = "ward.D2")
  order_receptors <- hclust_receptors$labels[hclust_receptors$order]

  dist_ligands <- stats::dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands <- stats::hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]

  order_receptors <- order_receptors %>%
    intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor <- order_ligands_receptor %>%
    intersect(colnames(lr_network_top_matrix))

  vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors,
                                                      order_ligands_receptor]
  if (send_spec == "mouse") {
    colnames(vis_ligand_receptor_network) <-
      colnames(vis_ligand_receptor_network) %>%
      nichenetr::convert_human_to_mouse_symbols() %>%
      as.character()
  }
  if (rec_spec == "mouse") {
    rownames(vis_ligand_receptor_network) <-
      rownames(vis_ligand_receptor_network) %>%
      nichenetr::convert_human_to_mouse_symbols() %>%
      as.character()
  }

  vis_ligand_receptor_network <-
    vis_ligand_receptor_network[!is.na(rownames(vis_ligand_receptor_network)),
                                !is.na(colnames(vis_ligand_receptor_network))]

  # If requested, create dot plots of the receptors and ligands
  if (d_plot == TRUE) {
    dotplot_reciever <- Seurat::DotPlot(subset(sobject, idents = receiver),
                                        features = rev(rownames(
                                          vis_ligand_receptor_network)),
                                        cols = c("gray95", "sienna3")) &
      Seurat::RotatedAxis()
    print(dotplot_reciever)

    dotplot_ligand <- Seurat::DotPlot(subset(sobject, idents = senders),
                                      features = rev(colnames(
                                        vis_ligand_receptor_network)),
                                      cols = c("gray95", "sienna3")) &
      Seurat::RotatedAxis() &
      ggplot2::coord_flip()
    print(dotplot_ligand)
  }

  # If requested, create a ligand-receptor visualization graph, weighted by
  # transcriptional changes
  if (lr_vis == TRUE) {
    lrv <- vis_ligand_receptor_network %>%
      t() %>%
      nichenetr::make_heatmap_ggplot("Ligands",
                          "Receptors",
                          x_axis_position = "top",
                          legend_title = "Prior interaction potential",
                          color = "darkseagreen3")
    print(lrv)
  } else {
    lrv <- NULL
  }

  return.list <- list(dotplot_ligand            = dotplot_ligand,
                      dotplot_receiver          = dotplot_reciever,
                      ligand_target_heatmap     = ltv,
                      ligand_receiver_heatmap   = lrv,
                      ligand_activities         = ligand_activities,
                      ligand_target_matrix      = vis_ligand_target,
                      ligand_receptor_matrix    = vis_ligand_receptor_network)

}


#' Load ligand receptor data
#'
#' @return None
#' @keywords internal
#'
load_lig_receptor_data <- function() {
  package_dir <- find.package("rrrSingleCellUtils")

  if (is.null(rrr_env$ligands)) {
    if (file.exists(paste(package_dir, "/LRT_Reference.RData", sep = ""))) {
      load(paste(package_dir, "/LRT_Reference.RData", sep = ""))
      rrr_env$ligands <- ligands
      rrr_env$ligands_bona_fide <- ligands_bona_fide
      rrr_env$receptors <- receptors
      rrr_env$receptors_bona_fide <- receptors_bona_fide
      rrr_env$ligand_target_matrix <- ligand_target_matrix
      rrr_env$lr_network <- lr_network
      rrr_env$lr_network_strict <- lr_network_strict
      rrr_env$weighted_networks <- weighted_networks
      rrr_env$weighted_networks_lr <- weighted_networks_lr
    } else {
      gen_lig_receptor_ref()
    }
  }
}

#' Get reference data for find_ligand() function
#'
#' @param lig_tar_matrix Ligand target matrix R data file location
#' @param lig_rec_network Ligand receptor network R data file location
#' @param weighted_network Ligand receptor weighted network R data file location
#'
#' @keywords internal
#' @return None
#'
gen_lig_receptor_ref <- function(
  lig_tar_matrix =
    "https://zenodo.org/record/3260758/files/ligand_target_matrix.rds",
  lig_rec_network =
    "https://zenodo.org/record/3260758/files/lr_network.rds",
  weighted_network =
    "https://zenodo.org/record/3260758/files/weighted_networks.rds"
  ) {

  message("Getting reference data for find_ligands(). This will download nearly
        1Gb of data")

  package_dir <- find.package("rrrSingleCellUtils")

  # Get raw data
  ligand_target_matrix <- readRDS(url(lig_tar_matrix))

  lr_network <- readRDS(url(lig_rec_network))

  ligands <- lr_network %>%
    dplyr::pull(from) %>%
    unique()

  receptors <- lr_network %>%
    dplyr::pull(to) %>%
    unique()

  lr_network_strict <- lr_network %>%
    dplyr::filter(database != "ppi_prediction_go" &
                    database != "ppi_prediction")

  ligands_bona_fide <- lr_network_strict %>%
    dplyr::pull(from) %>%
    unique()

  receptors_bona_fide <- lr_network_strict %>%
    dplyr::pull(to) %>%
    unique()

  weighted_networks <- readRDS(url(weighted_network))

  weighted_networks_lr <- weighted_networks$lr_sig %>%
    dplyr::inner_join(lr_network %>%
                        dplyr::distinct(from, to),
                      by = c("from", "to"))


  save(ligands, ligands_bona_fide, receptors, receptors_bona_fide,
       ligand_target_matrix, lr_network, lr_network_strict, weighted_networks,
       weighted_networks_lr,
       file = paste(package_dir, "/LRT_Reference.RData", sep = ""))

  # load the correct data into the rrr_env variables
  rrr_env$ligands <- ligands
  rrr_env$ligands_bona_fide <- ligands_bona_fide
  rrr_env$receptors <- receptors
  rrr_env$receptors_bona_fide <- receptors_bona_fide
  rrr_env$ligand_target_matrix <- ligand_target_matrix
  rrr_env$lr_network <- lr_network
  rrr_env$lr_network_strict <- lr_network_strict
  rrr_env$weighted_networks <- weighted_networks
  rrr_env$weighted_networks_lr <- weighted_networks_lr


  # Mouse data - #### Not currently used? ####
                 #### Seems like we convert mouse input to human genes ####
  # m_ligand_target_matrix <- ligand_target_matrix
  #
  # rownames(m_ligand_target_matrix) <-
  #   nichenetr::convert_human_to_mouse_symbols(rownames(m_ligand_target_matrix))
  # colnames(m_ligand_target_matrix) <-
  #   nichenetr::convert_human_to_mouse_symbols(colnames(m_ligand_target_matrix))
  #
  # m_ligand_target_matrix <-
  #   m_ligand_target_matrix[!is.na(rownames(m_ligand_target_matrix)),
  #                          !is.na(colnames(m_ligand_target_matrix))]
  #
  # m_lr_network <- lr_network
  # m_lr_network$from <-
  #   nichenetr::convert_human_to_mouse_symbols(m_lr_network$from)
  # m_lr_network$to <-
  #   nichenetr::convert_human_to_mouse_symbols(m_lr_network$to)
  # m_lr_network <- stats::na.omit(m_lr_network)
  #
  # m_ligands <- nichenetr::convert_human_to_mouse_symbols(ligands) %>%
  #   stats::na.omit() %>%
  #   as.character()
  #
  # m_receptors <- nichenetr::convert_human_to_mouse_symbols(receptors) %>%
  #   stats::na.omit() %>%
  #   as.character()
  #
  # m_ligands_bona_fide <-
  #   nichenetr::convert_human_to_mouse_symbols(ligands_bona_fide) %>%
  #   stats::na.omit() %>%
  #   as.character()
  #
  # m_receptors_bona_fide <-
  #   nichenetr::convert_human_to_mouse_symbols(receptors_bona_fide) %>%
  #   stats::na.omit() %>%
  #   as.character()
  #
  # m_lr_network_strict <- lr_network_strict
  #
  # m_lr_network_strict$from <-
  #   nichenetr::convert_human_to_mouse_symbols(m_lr_network_strict$from)
  # m_lr_network_strict$to <-
  #   nichenetr::convert_human_to_mouse_symbols(m_lr_network_strict$to)
  # m_lr_network_strict <- stats::na.omit(m_lr_network_strict)
  #
  # m_weighted_networks <- weighted_networks
  #
  # m_weighted_networks$lr_sig$from <-
  #   nichenetr::convert_human_to_mouse_symbols(m_weighted_networks$lr_sig$from)
  # m_weighted_networks$lr_sig$to <-
  #   nichenetr::convert_human_to_mouse_symbols(m_weighted_networks$lr_sig$to)
  # m_weighted_networks$lr_sig <- stats::na.omit(m_weighted_networks$lr_sig)
  #
  # m_weighted_networks$gr$from <-
  #   nichenetr::convert_human_to_mouse_symbols(m_weighted_networks$gr$from)
  # m_weighted_networks$gr$to <-
  #   nichenetr::convert_human_to_mouse_symbols(m_weighted_networks$gr$to)
  # m_weighted_networks$gr <- stats::na.omit(m_weighted_networks$gr)
  #
  # m_weighted_networks_lr <- m_weighted_networks$lr_sig %>%
  #   dplyr::inner_join(m_lr_network %>%
  #                       dplyr::distinct(from, to),
  #                     by = c("from", "to"))
  #
  # save(m_ligands, m_ligands_bona_fide, m_receptors, m_receptors_bona_fide,
  #      m_ligand_target_matrix, m_lr_network, m_lr_network_strict,
  #      m_weighted_networks, m_weighted_networks_lr,
  #      file = paste(package_dir, "/data/m_LRT_Reference.RData"))
}

# Make an environment with the variables inside to allow them to be modified
# by functions within the package down the line.
rrr_env <- new.env(parent = emptyenv())

rrr_env$ligands <- NULL
rrr_env$ligands_bona_fide <- NULL
rrr_env$receptors <- NULL
rrr_env$receptors_bona_fide <- NULL
rrr_env$ligand_target_matrix <- NULL
rrr_env$lr_network <- NULL
rrr_env$lr_network_strict <- NULL
rrr_env$weighted_networks <- NULL
rrr_env$weighted_networks_lr <- NULL
