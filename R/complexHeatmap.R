#' Make a complex heatmap showing nichenet output from FindLigands()
#'
#' @param fl_object An object returned by the FindLigands() function
#'
#' @return a ggplot object
#' @export
#'
#' @examples
plot_complex_heatmap <- function(fl_object) {

  # middle_plot <- fl_object[[4]]
  # right_plot <- fl_object[[1]]
  # bottom_plot <- fl_object[[2]]
  #

  mp_data <- fl_object[[4]]$data %>%
    mutate(x = factor(x, ordered = FALSE),
           y = factor(y, ordered = FALSE))

  rp_data <- fl_object[[1]]$data %>%
    as_tibble() %>%
    rename(x = id,
           y = features.plot)

  bp_data <- fl_object[[2]]$data %>%
    rename(x = features.plot,
           y = id)

  buffer_data <- tibble(x = factor(" "),
                        y = factor(" "))

  ##### Bind rows instead
  ###### put in a column for x-axis label colors above here and use
  ### axis.text.y = element_text(colour = c[[col]]) to plot
  combined <- bind_rows(mp_data, rp_data) %>%
    bind_rows(., bp_data) %>%
    bind_rows(buffer_data) %>%
    mutate(x = fct_relevel(x, c(levels(mp_data$x), levels(buffer_data$x), levels(rp_data$x)))) %>%
    mutate(y = fct_relevel(y, c(levels(mp_data$y), levels(buffer_data$y), levels(bp_data$y))))

  ggplot(combined, aes(x = x, y = y)) +
    geom_tile(aes(fill = score), width = 0.9, height = 0.9) +
    scale_fill_gradient(low = "gray95",
                        high = "midnightblue",
                        na.value="transparent",
                        name = "Int Score") +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_color_gradient(low = "white",
                         high = "darkred",
                         na.value="transparent",
                         name = "Scaled Exp") +
    scale_size(name = "% Exp") +
    xlab("Receptors") +
    ylab("Ligands") +
    guides(fill = guide_legend(order = 1),
           size = guide_legend(order = 2),
           shape = guide_legend(order = 3)) +
    theme(axis.text.x =
            element_text(angle = 90,
                         hjust = 0.9,
                         vjust = 0.6))
}
