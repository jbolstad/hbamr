#' Plot estimated respondent positions
#'
#' Plot the distribution of estimated respondent positions from an HBAM model.
#'
#' @export
#' @param object An instance of class `stanfit` produced by `hbam()`.
#' @param inc_stimuli Logical: Should estimated stimulus positions also be shown?
#' @param n_draws Integer specifying the number of posterior draws to use when illustrating the uncertainty of the population distribution. Defaults to 15.
#' @param color Color of lines illustrating uncertainty.
#' @param fill Fill color for density plots.
#' @param alpha_color Number in \[0,1\]: Inverse level of transparency for line color.
#' @param alpha_fill Number in \[0,1\]: Inverse level of transparency for fill color.
#' @param seed A positive integer specifying an optional seed for reproducibility. The seed is used to select respondent position draws for illustrating uncertainty.
#' @return A `ggplot` object.
#'

plot_respondents <- function(object, inc_stimuli = TRUE, n_draws = 15, color = "#053061", fill = "#2166AC", alpha_color = .6, alpha_fill = .7 / n_draws, seed = 1) {
  pd <- get_plot_data(object, n_draws = n_draws, seed = seed)
  suppressWarnings(
  if (inc_stimuli == T) {
    p <- ggplot2::ggplot() + ggplot2::geom_density(data = pd$chi_draws, ggplot2::aes(.data$chi, by = .data$draw_no), color = ggplot2::alpha(color, alpha_color), fill = ggplot2::alpha(fill, alpha_fill), na.rm = TRUE) +
      ggplot2::geom_text(data = pd$s_label, aes(x = .data$x, y = 0, label = .data$name), hjust = 0, nudge_y = .008, check_overlap = TRUE, angle = 90) +
      ggplot2::geom_point(data = pd$s_label, aes(x = .data$x, y = 0), color = color) +
      ggplot2::labs(x = "Ideological scale", y = "Density")
  } else {
    p <- ggplot2::ggplot() + ggplot2::geom_density(data = pd$chi_draws, ggplot2::aes(.data$chi, by = .data$draw_no), color = ggplot2::alpha(color, alpha_color), fill = ggplot2::alpha(fill, alpha_fill), na.rm = TRUE) +
      ggplot2::labs(x = "Ideological scale", y = "Density")
  })
  lim <- max(abs(quantile(pd$chi_draws$chi, probs = c(.002, .998))))
  p <- p + xlim(-lim, lim)
  return(p)
}
