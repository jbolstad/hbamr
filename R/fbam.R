#' Fit the FBAM_MINI model using optimization
#'
#' Fit a simplified Bayesian Aldrich-McKelvey model with fixed hyperparameters using optimization via `rstan`.
#'
#' @export
#' @param self A numerical vector of N ideological self-placements. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param stimuli An N × J matrix of numerical stimulus placements, where J is the number of stimuli. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param allow_miss Integer specifying how many missing stimulus positions to be accepted for an individual still to be included in the analysis. This argument will not be used if the data have been prepared in advance via the `prep_data` function. Defaults to 2.
#' @param req_valid Integer specifying how many valid observations to require for a respondent to be included in the analysis. The default is `req_valid = J - allow_miss`, but if specified, `req_valid` takes precedence. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param req_unique Integer specifying how may unique positions on the ideological scale each respondent is required to have used when placing the stimuli in order to be included in the analysis. The default is `req_unique = 2`. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param prep_data Logical: Should the data be prepared before fitting the model? (Or have the data been prepared in advance via the `prep_data` function? If so, set `prep_data = FALSE`.)
#' @param data List of data that have been prepared in advance via the `prep_data` function. Only applicable when `prep_data = TRUE`.
#' @param seed A positive integer specifying an optional seed for reproducibility. If this argument is not supplied, a random seed will be generated and the function will produce slightly different results on each run.
#' @param ... Arguments passed to `rstan::optimizing`.
#' @return A list produced by `rstan::optimizing`.
#'
#' @examples
#' \donttest{
#' # Loading ANES 2012 data:
#' data(LC2012)
#'
#' self <- LC2012[, 2]
#' stimuli <- LC2012[, -c(1:2)]
#'
#' # Fitting the FBAM_MINI model:
#' fit_fbam_mini <- fbam(self, stimuli)
#'
#' # Obtaining point estimates for the latent stimulus positions:
#' theta_est <- get_est(fit_fbam_mini, par = "theta")
#' }

fbam <- function(self = NULL, stimuli = NULL, allow_miss = 2, req_valid = NA,
                 req_unique = 2, prep_data = TRUE, data = NULL,
                 seed = sample.int(.Machine$integer.max, 1), ...) {
  if (prep_data == TRUE) { dat <- hbamr::prep_data(self, stimuli, allow_miss = allow_miss, req_valid = req_valid, req_unique = req_unique) } else { dat <- data }
  set.seed(seed)
  out <- rstan::optimizing(stanmodels[["FBAM_MINI"]], data = dat, init = inits[["FBAM_MINI"]](1, dat), seed = seed, ...)
  return(out)
}
