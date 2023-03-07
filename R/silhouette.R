#' Perform silhouette scoring on a Seurat object.
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

#' Find optimal FindClusters() resolution value to maximize silhouette score
#'
#' @details This function takes in a Seruat object returns a data frame of
#' silhouette scores for each resolution value
#' @param sobject A Seurat object containing all of the cells for analysis
#'    (required)
#' @param test_res A numeric vector of resolution values to test
#' @param summary_plot A logical value indicating whether to plot the results
#' @export
#'
#' @return A numeric value
#'
#' @examples
#' \dontrun{
#' optimize_silhouette(sobject = seurat_obj,
#'                     test_res = seq(0.1, 0.9, by = 0.1))
#' }
optimize_silhouette <- function(sobject,
                                test_res = seq(0.05, 0.75, by = 0.05),
                                summary_plot = TRUE) {

    if (.Platform$OS.type == "unix") {
        num_cores <- parallel::detectCores()
        output <-
            parallel::mclapply(test_res,
                               function(x) silhouette_to_df(sobject = sobject,
                                                            res = x),
                               mc.cores = num_cores) %>%
            dplyr::bind_rows()
    } else {
        output <-
            lapply(test_res, function(x) silhouette_to_df(sobject = sobject,
                                                          res = x)) %>%
            dplyr::bind_rows()
    }

    if (summary_plot) {
        print(ggplot2::ggplot(output,
                            ggplot2::aes(x = sil_vals,
                                         y = num_clusters)) +
            ggplot2::geom_text(ggplot2::aes(label = res_vals)) +
            ggplot2::labs(x = "Silhouette score",
                          y = "Number of clusters"))
    }

    return(output)
}

# Function to use with apply in optimize_silhouette instead of a for loop
silhouette_to_df <- function(sobject, res) {
    output <- list()

    output[["res_vals"]] <- res

    sobject <- sobject %>%
        Seurat::FindClusters(resolution = res, verbose = FALSE)
    if (max(as.numeric(as.character(sobject$seurat_clusters))) > 0) {
        sil_obj <- silhouette_seurat(sobject)

        output[["num_clusters"]] <- summary(sil_obj)$clus.avg.widths %>%
            length()
        output[["sil_vals"]] <- silhouette_mean(sil_obj)
    } else {
        output[["num_clusters"]] <- 1
        output[["sil_vals"]] <- NA
    }
    return(as.data.frame(output))
}
