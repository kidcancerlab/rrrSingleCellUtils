#' Make a complex heatmap showing nichenet output from FindLigands()
#'
#' @param fl_object An object returned by the FindLigands() function
#' @param grid_color_low Color for low values in the grid
#' @param grid_color_high Color for high values in the grid
#' @param point_color_low Color for low values in the points
#' @param point_color_high Color for high values in the points
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' data_obj <- findLigands(other_data)
#' plot_complex_heatmap(data_obj)
#' }
plot_complex_heatmap <- function(fl_object,
                                 grid_color_low = "gray95",
                                 grid_color_high = "midnightblue",
                                 point_color_low = "white",
                                 point_color_high = "darkred") {
  mp_data <- fl_object[[4]]$data %>%
    dplyr::mutate(x = factor(x, ordered = FALSE),
                  y = factor(y, ordered = FALSE))

  rp_data <- fl_object[[1]]$data %>%
    tibble::as_tibble() %>%
    dplyr::rename(x = id,
                  y = features.plot)

  bp_data <- fl_object[[2]]$data %>%
    dplyr::rename(x = features.plot,
                  y = id)

  buffer_data <- tibble::tibble(x = factor(" "),
                                y = factor(" "))

  ##### Bind rows instead
  ###### put in a column for x-axis label colors above here and use
  ### axis.text.y = element_text(colour = c[[col]]) to plot
  combined <- dplyr::bind_rows(mp_data, rp_data) %>%
    dplyr::bind_rows(., bp_data) %>%
    dplyr::bind_rows(buffer_data) %>%
    dplyr::mutate(x = forcats::fct_relevel(x,
                                           c(levels(mp_data$x),
                                             levels(buffer_data$x),
                                             levels(rp_data$x)))) %>%
    dplyr::mutate(y = forcats::fct_relevel(y,
                                           c(levels(mp_data$y),
                                             levels(buffer_data$y),
                                             levels(bp_data$y))))

  complex_plot <- ggplot2::ggplot(combined, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = score), width = 0.9, height = 0.9) +
    ggplot2::scale_fill_gradient(low = grid_color_low,
                                 high = grid_color_high,
                                 na.value = "transparent",
                                 name = "Int Score") +
    ggplot2::geom_point(ggplot2::aes(size = pct.exp, color = avg.exp.scaled)) +
    ggplot2::scale_color_gradient(low = point_color_low,
                                  high = point_color_high,
                                  na.value = "transparent",
                                  name = "Scaled Exp") +
    ggplot2::scale_size(name = "% Exp") +
    ggplot2::xlab("Receptors") +
    ggplot2::ylab("Ligands") +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1),
                    size = ggplot2::guide_legend(order = 2),
                    shape = ggplot2::guide_legend(order = 3)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x =
                     ggplot2::element_text(angle = 90,
                                           hjust = 0.9,
                                           vjust = 0.6))

  return(complex_plot)
}
