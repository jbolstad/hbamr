#' Fit an FBAM model using optimization
#'
#' Fit a simplified Bayesian Aldrich-McKelvey model with fixed hyperparameters using optimization via `rstan`. Users may replace the default priors by supplying their own values for the hyperparameters.
#'
#' @export
#' @param self A numerical vector of N ideological self-placements. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param stimuli An N × J matrix of numerical stimulus placements, where J is the number of stimuli. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param model Character: Name of the model to be used. Defaults to FBAM_MINI. The available options are the three models with "FBAM" in their name. See the documentation for the `hbam()` function for descriptions of the models.
#' @param allow_miss Integer specifying how many missing stimulus positions to be accepted for an individual still to be included in the analysis. This argument will not be used if the data have been prepared in advance via the `prep_data` function. Defaults to 2.
#' @param req_valid Integer specifying how many valid observations to require for a respondent to be included in the analysis. The default is `req_valid = J - allow_miss`, but if specified, `req_valid` takes precedence. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param req_unique Integer specifying how may unique positions on the ideological scale each respondent is required to have used when placing the stimuli in order to be included in the analysis. The default is `req_unique = 2`. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param data List of data that have been prepared in advance via the `prep_data` function. Not required if the arguments `self` and `stimuli` are provided.
#' @param group_id Integer vector of length N identifying which group each respondent belongs to. The supplied vector should range from 1 to the total number of groups in the data, and all integers between these numbers should be represented in the supplied data. These data are only required by models with "MULTI" in their name and will be ignored when fitting other models.
#' @param seed A positive integer specifying an optional seed for reproducibility. If this argument is not supplied, a random seed will be generated and the function will produce slightly different results on each run.
#' @param sigma_alpha A positive numeric value specifying the standard deviation of the prior on the shift parameters in the FBAM_MINI model, or the standard deviation of the parameters' deviation from the group-means in FBAM_MULTI models. (This argument will be ignored by HBAM models.) Defaults to B / 4, where B measures the length of the survey scale as the number of possible placements on one side of the center.
#' @param sigma_beta A positive numeric value specifying the standard deviation of the prior on the logged stretch parameters in the FBAM_MINI model, or the standard deviation of the logged parameters' deviation from the group-means in FBAM_MULTI models. (This argument will be ignored by HBAM models.) Defaults to .35.
#' @param sigma_mu_alpha A positive numeric value specifying the standard deviation of the prior on the group-means of the shift parameters in MULTI-type models. Defaults to B / 5.
#' @param sigma_mu_beta A positive numeric value specifying the standard deviation of the prior on the group-means of the logged stretch parameters in MULTI-type models. Defaults to .3.
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
                 req_unique = 2, group_id = NULL, data = NULL,
                 seed = sample.int(.Machine$integer.max, 1),
                 sigma_alpha = NULL, sigma_beta = .35,
                 sigma_mu_alpha = NULL, sigma_mu_beta = .3, ...) {
  if (!model %in% c("FBAM_MINI", "FBAM_MULTI", "FBAM_MULTI_NF")) { stop(paste(model, "is not a valid model choice for optimization.")) }
  if (!is.null(data) & (!is.null(self) | !is.null(stimuli))) { message("Note: When pre-prepared data are supplied, other data arguments will be ignored.") }
  if (is.null(data) & (is.null(self) | is.null(stimuli))) { message("Note: Required data not supplied.") }
  if (is.null(data)) { dat <- hbamr::prep_data(self, stimuli, allow_miss = allow_miss, req_valid = req_valid, req_unique = req_unique, group_id = group_id) } else { dat <- data }
  if (grepl("MULTI", model) & is.null(dat$gg)) { stop("No group_id supplied for MULTI-type model.") }
  if (!grepl("MULTI", model) & !is.null(dat$gg)) { message("Note: The supplied group_id will not be used as the chosen model is not a MULTI-type model.") }

  if (is.null(sigma_alpha)) { sigma_alpha <- dat$B / 4.0 }
  if (is.null(sigma_mu_alpha)) { sigma_mu_alpha <- dat$B / 5.0 }
  dat$sigma_alpha <- sigma_alpha
  dat$sigma_mu_alpha <- sigma_mu_alpha
  dat$sigma_beta <- sigma_beta
  dat$sigma_mu_beta <- sigma_mu_beta

  set.seed(seed)
  out <- rstan::optimizing(stanmodels[[model]], data = dat, init = inits[[model]](1, dat), seed = seed, ...)
  return(out)
}
