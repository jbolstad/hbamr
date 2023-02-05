#' Fit an HBAM model
#'
#' Fit a Hierarchical Bayesian Aldrich-McKelvey model using automatically tuned Hamiltonian Monte Carlo (NUTS) sampling via `rstan`.
#'
#' @export
#' @param self A numerical vector of N ideological self-placements. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param stimuli An N × J matrix of numerical stimulus placements, where J is the number of stimuli. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param model Character: Name of the model to be used. One of: `"HBAM"`, `"HBAM_2"`, `"HBAM_NE"`, `"HBAM_HM"`, `"HBAM_MINI"`, `"HBAM_0"`, `"HBAM_R"`, `"HBAM_R_MINI"`, or `"BAM"`. Defaults to `"HBAM"`.
#' @param allow_miss Integer specifying how many missing stimulus positions to be accepted for an individual still to be included in the analysis. This argument will not be used if the data have been prepared in advance via the `prep_data` function. Defaults to 2.
#' @param req_valid Integer specifying how many valid observations to require for a respondent to be included in the analysis. The default is `req_valid = J - allow_miss`, but if specified, `req_valid` takes precedence. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param req_unique Integer specifying how may unique positions on the ideological scale each respondent is required to have used when placing the stimuli in order to be included in the analysis. The default is `req_unique = 2`. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param prep_data Logical: Should the data be prepared before fitting the model? (Or have the data been prepared in advance via the `prep_data` function? If so, set `prep_data = FALSE`.)
#' @param data List of data that have been prepared in advance via the `prep_data` function. Only applicable when `prep_data = TRUE`.
#' @param prefs An N × J matrix of numerical stimulus ratings or preference scores. These data are only required by the `"HBAM_R"` and `"HBAM_R_MINI"` models and will be ignored when fitting other models.
#' @param chains A positive integer specifying the number of Markov chains. Defaults to 4.
#' @param cores The number of cores to use when executing the Markov chains in parallel. By default, all detected physical cores will be used if `chains` is equal to or higher than the number of cores.
#' @param warmup A positive integer specifying the number of warmup (aka burn-in) iterations per chain. If step-size adaptation is on (which it is by default), this also controls the number of iterations for which adaptation is run (and hence these warmup samples should not be used for inference). The number of warmup iterations should be smaller than `iter`.
#' @param iter A positive integer specifying the number of iterations for each chain (including warmup).
#' @param thin A positive integer specifying the period for saving samples.
#' @param control A named list of parameters to control the sampler's behavior. See the details in the documentation for the control argument in the `stan` function in the `rstan` package.
#' @param seed A positive integer specifying an optional seed for reproducibility. If this argument is not supplied, a random seed will be generated and the function will produce slightly different results on each run.
#' @param ... Arguments passed to `rstan::sampling`.
#' @details This package provides several alternative models, which can be specified using the names below. Users who are unsure which model to use are adviced to use the default HBAM model. If speed or sampling diagnostics are an issue, HBAM_MINI may provide a useful alternative.
#'
#' **HBAM** is the default model, which allows for scale flipping and employs hierarchical priors on all individual level parameters. It also models heteroskedastic errors that vary by both individual and stimuli. This model was introduced by Bølstad (2023).
#'
#' **HBAM_2** uses different hyperpriors for the shifting parameters of respondents with different self-placements. In particular, the model estimates a separate mean hyperparameter for each self-placement. This model avoids shrinking the shifting parameters toward a common population mean, and may therefore fit better than HBAM if there are clear differences in average shifting across the scale of self-placements. However, this model also tends to run slower than the standard model.
#'
#' **HBAM_NE** models the self-placements as if they contain no error. The latent respondent positions are not treated as parameters, but rather calculated as function of the self-placements and other individual level parameters. The respondents positions are not given a prior, which means the model relies on the likelihood function and the prior on beta to yield meaningful results for these positions. By assuming no error in the self-placements, the model may underestimate the uncertainty in estimated respondents positions, while otherwise yielding very similar results to the standard HBAM model. In contrast to the standard model, the estimated respondent positions from this model will not exhibit any shrinkage, which for some purposes may be desirable, as the results may better represent the true distances between respondents and stimuli. This model also runs somewhat faster than the standard HBAM model.
#'
#' **HBAM_HM** assumes the prediction errors in the stimuli placements to be homoskedastic. This simplified model should normally produce very similar results to the HBAM model, and it runs somewhat faster.
#'
#' **HBAM_MINI** combines the characteristics of HBAM_NE and HBAM_HM: It models the self-placements as if they contain no error and assumes the prediction errors in the stimuli placements to be homoskedastic. This is the simplest model provided in this package that still retains all key features of the HBAM model. It is also the fastest HBAM variant in this package -- sampling about twice as fast as the standard HBAM model for the dataset analyzed here (while yielding very similar point estimates). For large datasets, this model may provide a reasonable compromise between model complexity and estimation speed.
#'
#' **HBAM_0** does not allow for scale flipping. This may be useful if there are truly zero cases of scale flipping in the data. Such scenarios can be created artificially, but may also arise in real data. For example, expert surveys appear unlikely to contain many instances of scale flipping.
#'
#' **HBAM_R** incorporates the rationalization component of the ISR model by Bølstad (2020). This model requires additional data to be supplied to the `prep_data()` function: An N × J matrix of stimuli ratings from the respondents, supplied as the argument `pref`. The rationalization part of the model is simplified relative to the original ISR model: The direction in which respondents move disfavored stimuli is estimated as a common expectation for each possible self-placement on the scale.
#'
#' **HBAM_R_MINI** combines the features of the HBAM_R model with the light-weight features of the HBAM_MINI model to achieve faster sampling compared to HBAM_R.
#'
#' **BAM** is an unpooled model, similar to the JAGS version introduced by Hare et al. (2015). This model is mainly provided to offer a baseline for model comparisons. While it is simple and fast, this model tends to overfit the data and produce invalid posterior distributions for some respondent positions (Bølstad 2023).
#'
#' Some of these models can also be used in situations where self-placements are not available and the only goal is to estimate stimulus positions. This can be achieved by supplying a vector of zeros (or random data) instead of real self-placements: `self = rep(0, nrow(stimuli))`. The HBAM_NE and HBAM_MINI models are then the relevant alternatives, as the other HBAM variants will include superfluous parameters (and will not sample properly with zero variance in the supplied self-placement data).
#' @return An object of S4 class `stanfit`.

#' @references
#' - Bølstad, Jørgen (forthcoming). Hierarchical Bayesian Aldrich-McKelvey Scaling. <i>Political Analysis</i>.
#' - Bølstad, Jørgen (2020). Capturing Rationalization Bias and Differential Item Functioning: A Unified Bayesian Scaling Approach. <i>Political Analysis</i> 28(3): 340–355.
#' - Hare, Christopher et al. (2015). Using Bayesian Aldrich-McKelvey Scaling to Study Citizens' Ideological Preferences and Perceptions. <i>American Journal of Political Science</i> 59(3): 759–774.
#'
#' @examples
#' # Loading and re-coding ANES 1980 data:
#' data(LC1980)
#' dat <- LC1980
#' dat[dat == 0 | dat == 8 | dat == 9] <- NA
#'
#' # Making a small subset of the data for illustration:
#' self <- dat[1:50, 1]
#' stimuli <- dat[1:50, -1]
#'
#' # Fitting the HBAM_MINI model:
#' fit_hbam <- hbam(self, stimuli, model = "HBAM_MINI",
#'                  warmup = 300, iter = 500, chains = 2, thin = 1)
#'
#' # Plot distribution of estimated respondent positions:
#' plot_respondents(fit_hbam)
#'
#' \dontrun{
#' # Fitting the standard HBAM model to the complete ANES 1980 data:
#' self <- dat[, 1]
#' stimuli <- dat[, -1]
#' fit_hbam <- hbam(self, stimuli)
#'
#' # Preparing the data before fitting, using defaults:
#' dat <- prep_data(self, stimuli)
#' fit_hbam <- hbam(data = dat, prep_data = FALSE)
#' }

hbam <- function(self = NULL, stimuli = NULL, model = "HBAM", allow_miss = 2, req_valid = NA,
                 req_unique = 2, prefs = NULL, prep_data = TRUE, data = NULL,
                 chains = 4, cores = parallel::detectCores(logical = FALSE),
                 warmup = 1000, iter = 4000, thin = 3, control = list(adapt_delta = .6),
                 seed = sample.int(.Machine$integer.max, 1), ...) {
  if (prep_data == TRUE) { dat <- hbamr::prep_data(self, stimuli, prefs, allow_miss = allow_miss, req_valid = req_valid, req_unique = req_unique) } else { dat <- data }
  set.seed(seed)
  init_ll <- lapply(1:chains, function(id) inits[[model]](id, dat))
  out <- rstan::sampling(stanmodels[[model]], data = dat, init = init_ll,
                         chains = chains, cores = cores, warmup = warmup, iter = iter, thin = thin, control = control, seed = seed, ...)
  return(out)
}
