#' Extract summaries of marginal posterior distributions from an HBAM model or point estimates from an FBAM model
#'
#' For objects produced by `hbam`, this function is a wrapper for `rstan::summary`. For objects produced by `fbam` it offers a way to extract point estimates.
#'
#' @export
#' @param object An instance of class `stanfit` produced by `hbam`, or a list produced by `fbam`.
#' @param par Character: Name of the parameter type to be extracted. Typically `"theta"` (stimuli positions) or `"chi"` (respondent positions).
#' @param probs A numeric vector of quantiles of interest for summarizing `stanfit` objects. The default is c(0.025, 0.50, 0.975).
#' @param simplify Logical: Should the returned object be simplified by dropping the Monte Carlo standard error and the posterior standard deviation? Defaults to `TRUE`.
#' @param ... Other arguments are passed on to `rstan::summary` when summarizing `stanfit` objects.
#' @return A tibble containing summaries of marginal posterior distributions. For objects produced by `fbam`, only maximum a posteriori estimates are returned.

get_est <- function (object, par = "theta", probs = c(0.025, 0.50, 0.975), simplify = TRUE, ...) {
  if (inherits(object, "stanfit")) {
    out <- dplyr::as_tibble(rstan::summary(object, par, probs = probs, ...)[[1]])
    if (simplify == TRUE) { out <- out[, -c(2, 3)] }
  } else {
    if (inherits(object, "list")) {
      out <- dplyr::as_tibble(object$par[startsWith( names(object$par), paste0(par, "[") )])
      names(out) <- par
    }
  }

  return(out)
}
