#' Create a gene list containing putative targets of ligand activity.
#'
#' @details This function takes in a Seruat object and runs silhouette scoring
#' @param sobject A Seurat object containing all of the cells for analysis
#'    (required)
#' @param cluster_col A column name containing the cluster assignments for cells
#' @param dims Numeric vector of PCA dimensions to use
#' @export
#'
#' @return A silhouette object
#'
#' @examples
#' \dontrun{
#' sil <- silhouette_seurat(sobject = seurat_obj)
#' }
silhouette_seurat <- function(sobject,
                              cluster_col = "seurat_clusters",
                              dims = c(1:10)) {
    sil_obj <- cluster::silhouette(
        x = sobject@meta.data[[cluster_col]] %>%
            as.character() %>%
            as.numeric(),
        dist = Seurat::Embeddings(sobject,
                                  reduction = "pca")[, dims] %>%
        stats::dist())

    return(sil_obj)
}

#' Calculate the average silhouette score from a silhouette object
#'
#' @details This function takes in a Seruat object and runs silhouette scoring
#' @param sobject A Seurat object containing all of the cells for analysis
#'    (required)
#' @param sil_obj A silhouette object
#' @export
#'
#' @return A numeric value
#'
#' @examples
#' \dontrun{
#' sil <- silhouette_seurat(sobject = seurat_obj)
#' sil_mean <- silhouette_mean(sil_obj = sil)
#' }
silhouette_mean <- function(sil_obj) {
    sil_summary <- summary(sil_obj)
    return(sil_summary$avg.width)
}

#' Calculate the average silhouette score from each cluster in a silhouette object
#'
#' @details This function takes in a silhouette object and returns a vector of
#' mean silhouette scores for each cluster
#' @param sobject A Seurat object containing all of the cells for analysis
#'    (required)
#' @param sil_obj A silhouette object
#' @export
#'
#' @return A numeric value
#'
#' @examples
#' \dontrun{
#' sil <- silhouette_seurat(sobject = seurat_obj)
#' sil_cluster_mean <- silhouette_cluster_mean(sil_obj = sil)
#' }
silhouette_cluster_mean <- function(sil_obj) {
    sil_summary <- summary(sil_obj)
    return(sil_summary$clus.avg.widths)
}

#' Plot silhouette scores
#'
#' @details This function takes in a Seruat object and runs silhouette scoring
#' @param sobject A Seurat object containing all of the cells for analysis
#'    (required)
#' @param sil_obj A silhouette object
#' @export
#'
#' @return A numeric value
#'
#' @examples
#' \dontrun{
#' sil <- silhouette_seurat(sobject = seurat_obj)
#' silhouette_plot(sil_obj = sil)
#' }
silhouette_plot <- function(sil_obj) {
    widths <- sil_obj[, ] %>%
        as.data.frame() %>%
        dplyr::group_by(cluster) %>%
        dplyr::arrange(sil_width, .by_group = TRUE) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(order = 1:nrow(sil_obj))

    plot_name <- ggplot2::ggplot(widths,
                                 ggplot2::aes(x = order,
                                              y = sil_width,
                                              fill = as.factor(cluster))) +
     ggplot2::geom_col() +
     ggplot2::labs(x = "", y = "Silhouette width") +
     ggplot2::scale_fill_discrete(name = "Cluster") +
     ggplot2::coord_flip()

    return(plot_name)
}
