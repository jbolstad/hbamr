#' Fit an FBAM model using optimization
#'
#' Fit a simplified Bayesian Aldrich-McKelvey model with fixed hyperparameters using optimization via `rstan`.
#'
#' @export
#' @param self A numerical vector of N ideological self-placements. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param stimuli An N × J matrix of numerical stimulus placements, where J is the number of stimuli. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param model Character: Name of the model to be used. Defaults to FBAM_MINI. The options are the three models with "FBAM" in the name. See the documentation for the `hbam()` function for descriptions of the models.
#' @param allow_miss Integer specifying how many missing stimulus positions to be accepted for an individual still to be included in the analysis. This argument will not be used if the data have been prepared in advance via the `prep_data` function. Defaults to 2.
#' @param req_valid Integer specifying how many valid observations to require for a respondent to be included in the analysis. The default is `req_valid = J - allow_miss`, but if specified, `req_valid` takes precedence. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param req_unique Integer specifying how may unique positions on the ideological scale each respondent is required to have used when placing the stimuli in order to be included in the analysis. The default is `req_unique = 2`. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param prep_data Logical: Should the data be prepared before fitting the model? (Or have the data been prepared in advance via the `prep_data` function? If so, set `prep_data = FALSE`.)
#' @param data List of data that have been prepared in advance via the `prep_data` function. Only applicable when `prep_data = TRUE`.
#' @param group_id Integer vector of length N identifying which group each respondent belongs to. The supplied vector should range from 1 to the total number of groups in the data, and all integers between these numbers should be represented in the supplied data. These data are only required by models with "MULTI" in their name and will be ignored when fitting other models.
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

fbam <- function(self = NULL, stimuli = NULL, model = "FBAM_MINI", allow_miss = 2, req_valid = NA,
                 req_unique = 2, group_id = NULL, prep_data = TRUE, data = NULL,
                 seed = sample.int(.Machine$integer.max, 1), ...) {
  if (!model %in% c("FBAM_MINI", "FBAM_MULTI", "FBAM_MULTI_0")) {
    stop(paste(model, "is not a valid model choice for optimization."))
  }
  if (prep_data == TRUE) { dat <- hbamr::prep_data(self, stimuli, allow_miss = allow_miss, req_valid = req_valid, req_unique = req_unique, group_id = group_id) } else { dat <- data }
  if (grepl("MULTI", model) & is.null(dat$gg)) { stop("No group_id supplied for MULTI-type model.") }
  set.seed(seed)
  out <- rstan::optimizing(stanmodels[[model]], data = dat, init = inits[[model]](1, dat), seed = seed, ...)
  return(out)
}
