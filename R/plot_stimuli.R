#' Plot estimated stimulus positions
#'
#' Plot marginal posterior distributions of stimulus positions from an HBAM model
#'
#' @export
#' @param object An instance of class `stanfit` produced by `hbam`.
#' @param rev_color Logical: Display low positions as red and high positions as blue.
#' @param alpha Number in \[0,1\]: Inverse level of transparency for fill color.
#' @return A `ggplot` object.
#'

plot_stimuli <- function(object, rev_color = FALSE, alpha = .55) {
  pd <- get_plot_data(object)
  pal <- RColorBrewer::brewer.pal(n = 11, name = "RdBu")
  if (rev_color == F) { pal <- rev(pal) }
  # Placing stimuli in 11 equally spaced categories for fill-color:
  cuts = (-5:6 - .5) * max(abs(pd$s_label$x)) * 2 / 10
  cats <- as.numeric(cut(pd$s_label$x, breaks = cuts))

  p <- ggplot2::ggplot(pd$s_draws, ggplot2::aes(.data$value, fill = .data$name)) + ggplot2::geom_density(color = "gray50", alpha = alpha) +
    ggplot2::geom_text(data = pd$s_label, ggplot2::aes(x = .data$x, y = 0, label = .data$name), hjust = 0, nudge_x = -.007,
              nudge_y = .025 * max(pd$s_label$mode_y), check_overlap = TRUE, angle = 90) +
    ggplot2::labs(x = "Ideological scale", y = "Posterior density") +
    ggplot2::theme(legend.position = "none") + ggplot2::scale_fill_manual(values = pal[cats])
  return(p)
}

