#' Extract data for plotting results from an HBAM model
#'
#' @export
#' @param object An instance of class `stanfit`.
#' @param n_draws Integer specifying the number of posterior draws to use when illustrating the uncertainty of the population distribution.
#' @param seed A positive integer specifying an optional seed for reproducibility. The seed is used to select respondent position draws for illustrating uncertainty.
#' @return A list of three `tibble`s: The first element contains the posterior mean stimulus positions, as well as the x- and y-values of the posterior modes (which can be useful for labeling the distributions). The second element contains the posterior draws for the stimulus positions (which can be used to calculate marginal posterior densities). The third element contains the selected number of posterior draws for each respondent (which form the key ingredient for `plot_respondents`).
#'

get_plot_data <- function(object, n_draws = 25, seed = 1) {
  set.seed(seed)
  chi_draws <- as.data.frame((rstan::extract(object, pars = "chi")$chi)[1:n_draws, ])
  chi_draws$draw_no <- rownames(chi_draws)
  chi_draws <- pivot_longer(chi_draws, starts_with("V"), values_to = "chi")
  s_draws <- as.data.frame(rstan::extract(object, pars = "theta")$theta)

  # Use stimuli-names if available:
  if (!is.null(names(attr(object@sim$samples[[1]], "args")$init_list$theta_raw))) {
    colnames(s_draws) <- names(attr(object@sim$samples[[1]], "args")$init_list$theta_raw) }
  s_draws <- pivot_longer(s_draws, everything())

  label <- s_draws %>%
    dplyr::group_by(.data$name) %>%
    dplyr::summarize(x = mean(.data$value), mode_x = get_post_mode_x(.data$value), mode_y = get_post_mode_y(.data$value))
  # label$y_adj <- label$y + mean(label$y) * .07
  return(list(s_label = label, s_draws = s_draws, chi_draws = chi_draws))
}
