#' Extract data for plotting results from an HBAM model
#'
#' Extract data for plotting results from an HBAM model.
#'
#' @export
#' @param object An instance of class `stanfit` produced by `hbam()` or a list produced by `fbam()`.
#' @param n_draws Integer specifying the number of posterior draws to use when illustrating the uncertainty of the population distribution. This only applies for `stanfit` objects.
#' @param seed A positive integer specifying an optional seed for reproducibility. The seed is used to select respondent position draws for illustrating uncertainty. This only applies for `stanfit` objects.
#' @return A list of three `tibble`s: The first element contains the posterior mean stimulus positions, as well as the x- and y-values of the posterior modes (which can be useful for labeling the distributions). The second element contains the posterior draws for the stimulus positions (which can be used to calculate marginal posterior densities). The third element contains the selected number of posterior draws for each respondent (which form the key ingredient for `plot_respondents`).
#'

get_plot_data <- function(object, n_draws = 15, seed = 1) {
  if (inherits(object, "stanfit")) {
    set.seed(seed)
    chi_draws <- as.data.frame((rstan::extract(object, pars = "chi")$chi))
    chi_draws <- chi_draws[sample(nrow(chi_draws), n_draws), ]
    chi_draws$draw_no <- rownames(chi_draws)
    chi_draws <- pivot_longer(chi_draws, starts_with("V"), values_to = "chi")
    s_draws <- as.data.frame(rstan::extract(object, pars = "theta")$theta)

    # Use stimuli-names if available:
    if (!is.null(names(attr(object@sim$samples[[1]], "args")$init_list$theta_raw))) {
      colnames(s_draws) <- names(attr(object@sim$samples[[1]], "args")$init_list$theta_raw) }
    s_draws <- pivot_longer(s_draws, everything())
    s_draws$name <- factor(s_draws$name)

    label <- s_draws %>%
      dplyr::group_by(.data$name) %>%
      dplyr::summarize(x = mean(.data$value), mode_x = get_post_mode_x(.data$value), mode_y = get_post_mode_y(.data$value))
    # label$y_adj <- label$y + mean(label$y) * .07
  } else {
    if (inherits(object, "list")) {
      chi_draws <- dplyr::as_tibble(data.frame(chi = object$par[grepl( "chi[", names(object$par), fixed = TRUE)]))
      chi_draws$draw_no <- 1
      s_draws <- data.frame(theta = object$par[grepl( "theta[", names(object$par), fixed = TRUE)])
      label <- dplyr::as_tibble(data.frame(x = s_draws$theta, mode_x = s_draws$theta, mode_y = NA))
      label$name <- rownames(s_draws)
      s_draws <- dplyr::as_tibble(s_draws)
    }
  }
  return(list(s_label = label, s_draws = s_draws, chi_draws = chi_draws))
}
