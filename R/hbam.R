#' Fit an HBAM model
#'
#' Fit a Hierarchical Bayesian Aldrich-McKelvey model using automatically tuned Hamiltonian Monte Carlo (NUTS) sampling via `rstan`.
#'
#' @export
#' @param self A numerical vector of N ideological self-placements. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param stimuli An N × J matrix of numerical stimulus placements, where J is the number of stimuli. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param model Character: Name of the model to be used. One of: `"HBAM"`, `"HBAM_2"`, `"HBAM_0"`, `"HBAM_R"`, or `"BAM"`. Defaults to `"HBAM"`.
#' @param allow_miss Integer specifying how many missing stimulus positions to be accepted for an individual still to be included in the analysis. This argument will not be used if the data have been prepared in advance via the `prep_data` function. Defaults to 2.
#' @param req_valid Integer specifying how many valid observations to require for a respondent to be included in the analysis. The default is `req_valid = J - allow_miss`, but if specified, `req_valid` takes precedence. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param req_unique Integer specifying how may unique positions on the ideological scale each respondent is required to have used when placing the stimuli in order to be included in the analysis. The default is `req_unique = 2`. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param prep_data Logical: Should the data be prepared before fitting the model? (Or have the data been prepared in advance via the `prep_data` function? If so, set `prep_data = FALSE`.)
#' @param data List of data that have been prepared in advance via the `prep_data` function. Only applicable when `prep_data = TRUE`.
#' @param chains A positive integer specifying the number of Markov chains. Defaults to 4.
#' @param cores The number of cores to use when executing the Markov chains in parallel. Defaults to `min(parallel::detectCores(logical = FALSE), chains, 6)`.
#' @param warmup A positive integer specifying the number of warmup (aka burn-in) iterations per chain. If step-size adaptation is on (which it is by default), this also controls the number of iterations for which adaptation is run (and hence these warmup samples should not be used for inference). The number of warmup iterations should be smaller than `iter`.
#' @param iter A positive integer specifying the number of iterations for each chain (including warmup).
#' @param thin A positive integer specifying the period for saving samples.
#' @param ... Arguments passed to `rstan::sampling`.
#' @return An object of S4 class `stanfit`.
#'
#' @examples
#'
#' # Loading and recoding data
#' data(LC1980)
#' dat <- LC1980
#' dat[dat == 0 | dat == 8 | dat == 9] <- NA
#' self <- dat[, 1]
#' stimuli <- dat[, -1]
#'
#' # Fitting the standard HBAM model, using default settings:
#' fit_hbam <- hbam(self, stimuli)
#'
#' # Preparing data in advance, using defaults:
#' dat <- prep_data(self, stimuli)
#' fit_hbam <- hbam(data = dat, prep_data = FALSE)

hbam <- function(self = NULL, stimuli = NULL, model = "HBAM", allow_miss = 2, req_valid = NA,
                 req_unique = 2, prefs = NULL, prep_data = TRUE, data = NULL,
                 chains = 4, cores = min(parallel::detectCores(logical = FALSE), chains, 6),
                 warmup = 1000, iter = 4000, thin = 3, control = list(adapt_delta = .6), ...) {
  if (prep_data == TRUE) { dat <- hbamr::prep_data(self, stimuli, prefs, allow_miss = allow_miss, req_valid = req_valid, req_unique = req_unique) } else { dat <- data }
  init_ll <- lapply(1:chains, function(id) hbamr:::inits[[model]](id, dat))
  out <- rstan::sampling(hbamr:::stanmodels[[model]], data = dat, init = init_ll,
                         chains = chains, cores = cores, warmup = warmup, iter = iter, thin = thin, control = control, ...)
  return(out)
}
