#' Extract point estimates or other summaries of marginal posterior distributions
#'
#' For objects produced by `hbam()`, this function is a wrapper for `rstan::summary()`. For objects produced by `fbam()` it offers a way to extract point estimates.
#'
#' @export
#' @param object An instance of class `stanfit` produced by `hbam()`, or a list produced by `fbam()`.
#' @param par Character: Name of the parameter type to be extracted. Typically `"theta"` (stimuli positions) or `"chi"` (respondent positions).
#' @param format_orig Logical: Should individual-level parameters be mapped to the original dataset by returning rows of NAs for respondents who were not included in the analysis? Defaults to `FALSE`.
#' @param probs A numeric vector of quantiles of interest for summarizing `stanfit` objects.
#' @param simplify Logical: Should the returned object be simplified by dropping the Monte Carlo standard error and the posterior standard deviation? Defaults to `TRUE`.
#' @param ... Other arguments are passed on to `rstan::summary()` when summarizing `stanfit` objects.
#' @return A tibble containing summaries of marginal posterior distributions. For objects produced by `fbam()`, only maximum a posteriori estimates are returned.

get_est <- function (object, par = "theta", format_orig = FALSE, probs = c(0.025, 0.50, 0.975), simplify = TRUE, ...) {
  if (inherits(object, "stanfit")) {
    out <- as.data.frame(rstan::summary(object, par, probs = probs, ...)[[1]])
    if (simplify == TRUE) { out <- out[, -c(2, 3)] }
    data <- object@.MISC$hbam_data
  } else {
    if (inherits(object, "list")) {
      if (par %in% names(object$par)) {
        out <- matrix(object$par[which(names(object$par) == par)], ncol = 1)
      } else {
        out <- matrix(object$par[startsWith( names(object$par), paste0(par, "[") )], ncol = 1)
      }
      colnames(out) <- par
      data <- object$hbam_data
    }
  }
  if (format_orig) {
    if (is.null(data) || !inherits(data, "hbam_data")) {
      stop("Could not find the data used for fitting within the supplied object. You need to supply an object produced by hbam() or fbam().")
    }
    if (nrow(out) == data$N) {
      out_long <- matrix(NA, ncol = ncol(out), nrow = length(data$keep))
      out_long[data$keep, ] <- unlist(out[, ])
      colnames(out_long) <- colnames(out)
      out <- out_long
    }
  }
  out <- dplyr::as_tibble(out)

  return(out)
}
