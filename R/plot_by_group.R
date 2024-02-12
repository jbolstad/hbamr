#' Plot posterior densities of parameter averages by group
#'
#' Plot posterior densities of group summaries of individual parameters. The respondents can be grouped by any categorical variable and the function works whether the fitted model is of "MULTI"-type or not.
#'
#' @export
#' @param object An instance of class `stanfit` produced by `hbam()`.
#' @param data The list of data that was used to produce the object.
#' @param group_id A vector that will be used to split the respondents into groups before plotting. The vector must either be as long as the number of rows in the original dataset, or as long as the number of respondents included in the analysis. The `group_id` supplied here can be the same as the `group_id` argument used when fitting a "MULTI"-type model, but it does not have to be.
#' @param par Character: Name of the parameter to be plotted. One of the following: `"alpha"`, `"beta"`, `"abs_beta"`, `"lambda"`, or `"chi"`. Defaults to `"abs_beta"`, which means the absolute values of the draws for beta will be used. Further individual-level parameters like `"eta"` can be specified if these have been passed to `hbam()` via the argument `extra_pars` when fitting the model. (Note that homoskedastic models have no `"eta"` parameters and "NF"-type models have no `"lambda"` or `"kappa"` parameters.)
#' @param fill Fill color. Passed on to `ggplot2::geom_density()`.
#' @param color Color of outer lines. Passed on to `ggplot2::geom_density()`.
#' @param alpha Number in \[0,1\]: Inverse level of transparency.
#' @param ncol Number of columns. The default aims to have approximately ten subplots per column.
#' @param ascending Logical: Should the groups be placed in ascending order based on their posterior means (`TRUE`) or should they be ordered based on their names (`FALSE`)? Defaults to `TRUE`.
#' @return A `ggplot` object.
#'


plot_by_group <- function(object, data, group_id, par = "abs_beta", fill = "#2166AC", color = "#053061", alpha = .5, ncol = max(1, round(length(unique(group_id))/10)), ascending_means = TRUE) {
  if (!(length(group_id) == data$N_orig | length(group_id) == dat$N)) {
    stop("The supplied group_id must either be the same length as the original data or the number of respondents included in the analysis.") }
  if (length(group_id) == data$N_orig) {
    group_id <- group_id[dat$keep]
  }
  if (par == "abs_beta") {
    par <- "beta"
    draws <- abs(rstan::extract(object, pars = par)[[par]])
  } else {
    draws <- rstan::extract(object, pars = par)[[par]]
  }
  means <- apply(draws, 1, function(x) tapply(x, group_id, mean))
  mean_means <- apply(means, 1, mean)
  group_names <- attributes(means)$dimnames[[1]]
  means <- data.frame(group_id = group_names, means)
  if (ascending_means) {
    means$group_id <- factor(means$group_id, levels = group_names[order(mean_means)])
  } else {
    means$group_id <- factor(means$group_id, levels = group_names)
  }
  means_long <- tidyr::pivot_longer(means, !group_id, names_to = "draw_no")

  lim <- quantile(means_long$value, probs = c(.001, .999))
  plot <- ggplot2::ggplot(means_long, aes(x = value)) + ggplot2::geom_density(fill = fill, color = color, alpha = alpha, na.rm = TRUE) +
    ggplot2::xlim(lim[1], lim[2]) + ggplot2::labs(x = par, y = "Posterior density") + ggplot2::facet_wrap(~group_id, ncol = ncol)
  return(plot)
}

