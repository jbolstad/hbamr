#' Extract summaries of marginal posterior distributions from an HBAM model
#'
#' This function is a wrapper for wrapper for `rstan::summary`.
#'
#' @export
#' @param object An instance of class `stanfit` produced by `hbam`.
#' @param par Character: Name of the parameter type to be extracted. Typically `"theta"` (stimuli positions) or `"chi"` (respondent positions).
#' @param probs A numeric vector of quantiles of interest. The default is c(0.025, 0.50, 0.975).
#' @param simplify Logical: Should the returned object be simplified by dropping the Monte Carlo standard error and the posterior standard deviation? Defaults to `TRUE`.
#' @param ... Other arguments are passed on to `rstan::summary`.
#' @return A tibble containing summaries of marginal posterior distributions.

get_est <- function (object, par = "theta", probs = c(0.025, 0.50, 0.975), simplify = TRUE, ...) {
  out <- dplyr::as_tibble(rstan::summary(object, par, probs = probs, ...)[[1]])
  if (simplify == TRUE) { out <- out[, -c(2, 3)] }
  return(out)
}
