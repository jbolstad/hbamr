#' Plot individual parameter estimates over self-placements
#'
#' Create a boxplot of individual parameter estimates from an HBAM model over self-placements
#'
#' @export
#' @param object An object of class `stanfit` produced by `hbam`, or a `list` of such objects, which will produce a faceted plot.
#' @param data The list of data that was used to produce the `stanfit` object(s).
#' @param par Character: Name of the parameter to be plotted. One of the following: `"alpha"`, `"beta"`, `"abs_beta"`, `"lambda"`, `"chi"`, and `"eta"`. Defaults to `"chi"`.
#' @param estimate Character: Specifying which type of posterior point estimate to use. One of `"median"` and `"mean"`. Defaults to `"median"`.
#' @param names An optional character vector of model names of same length as the supplied list of models.
#' @param parlabel An optional character containing an alternative label for the parameter (will be parsed if passed as an expression).
#' @param fill Fill color of boxes. Passed on to `geom_boxplot`.
#' @param color Color of outer lines. Passed on to `geom_boxplot`.
#' @param width Width of boxes. Passed on to `geom_boxplot`.
#' @param alpha Number in \[0,1\]: Inverse level of transparency for fill color. Passed on to `geom_boxplot`.
#' @param outlier.size Size of dots representing outliers. Passed on to `geom_boxplot`.
#' @param median_color Color of solid line representing the median.
#' @param median_lwd Thickness of solid line representing the median.
#' @return A `ggplot` object.
#'

plot_over_self <- function(object, data, par = "chi", estimate = "median", names = NULL, parlabel = NULL,
                           fill = "#2166AC", color = "#053061", width = .7, alpha = .5, outlier.size = 0.3,
                           median_color = "black", median_lwd = .7) {
  if(is.null(parlabel)) { parlabel <- par}
  if(length(object) == 1) {
    pd <- get_pd(object, data, par, estimate)
    p <- ggplot2::ggplot(pd, ggplot2::aes(.data$V, .data$parameter)) + ggplot2::geom_boxplot(fill = fill, color = color, width = width, alpha = alpha, outlier.size = outlier.size) +
      xlab("Self-placement") + ylab(par)
    md <- ggplot2::ggplot_build(p)$data[[1]]
    p <- p + ggplot2::geom_segment(data = md, aes(x = .data$xmin, xend = .data$xmax,
                                   y = .data$middle, yend = .data$middle), colour = median_color, linewidth = median_lwd)
  } else {
    pd <- vector()
    for (m in 1:length(object)) {
      if (is.null(names)) { name <- object[[m]]@model_name } else { name <- names[m] }
      pd <- rbind(pd, dplyr::bind_cols(get_pd(object[[m]], data, par, estimate), model = rep(name, data$N)))
    }
      if (is.null(names)) {
        pd$model <- factor(pd$model, levels = unique(pd$model), labels = unique(pd$model)) } else {
        pd$model <- factor(pd$model, levels = names, labels = names) }
      p <- ggplot2::ggplot(pd, ggplot2::aes(.data$V, .data$parameter)) + ggplot2::geom_boxplot(fill = fill, color = color, width = width, alpha = alpha, outlier.size = outlier.size) +
        ggplot2::xlab("Self-placement") + ggplot2::ylab(par) +
        ggplot2::facet_wrap(~model, scales = "free") + ggplot2::ylab(parlabel)
      md <- ggplot2::ggplot_build(p)$data[[1]]
      md$model <- factor(md$PANEL, levels = as.numeric(unique(pd$model)), labels = levels(pd$model))
      p <- p + ggplot2::geom_segment(data = md, aes(x = .data$xmin, xend = .data$xmax,
                                            y = .data$middle, yend = .data$middle), colour = median_color, linewidth = median_lwd)
  }
  return(p)
}

get_pd <- function(object, data, par, estimate){
  if (par == "abs_beta") {
    draws <- as.matrix(rstan::extract(object, pars = "beta")$beta)
    param <- data.frame(
      median = apply(draws, 2, function(x) median(abs(x)) ),
      mean = apply(draws, 2, function(x) mean(abs(x)) ))
    names(param) <- c("50%", "mean")
  } else {
    param <- get_est(object, par)
    if (par == "eta") {
      param[, 1] <- sqrt(param[, 1]) / data$J # To get average sigma for each i over j.
      param[, 3] <- sqrt(param[, 3]) / data$J
      par <- "sqrt(eta) / J"
    }
  }
  if (estimate == "median") { pd <- dplyr::bind_cols(parameter = param$`50%`, V = as.ordered(data$V)) }
  if (estimate == "mean") { pd <- dplyr::bind_cols(parameter = param$mean, V = as.ordered(data$V)) }
  return(pd)
}
