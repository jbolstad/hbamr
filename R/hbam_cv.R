#' Perform K-fold cross-validation
#'
#' This function performs K-fold cross-validation for an HBAM or FBAM model in order to estimate the expected log pointwise predictive density for a new dataset (ELPD). Multiple chains for one or more folds can be run in parallel using the `future` package.
#'
#' @export
#' @param self A numerical vector of N ideological self-placements. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data()` function.
#' @param stimuli An N × J matrix of numerical stimulus placements, where J is the number of stimuli. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data()` function.
#' @param model Character: Name of the model to be used. Defaults to HBAM.
#' @param allow_miss Integer specifying how many missing stimulus positions should be accepted for an individual still to be included in the analysis. This argument will not be used if the data have been prepared in advance via the `prep_data()` function. Defaults to 0.
#' @param req_valid Integer specifying how many valid observations to require for a respondent to be included in the analysis. The default is `req_valid = J - allow_miss`, but if specified, `req_valid` takes precedence. This argument will not be used if the data have been prepared in advance via the `prep_data()` function.
#' @param req_unique Integer specifying how may unique positions on the ideological scale each respondent is required to have used when placing the stimuli in order to be included in the analysis. The default is `req_unique = 2`. This argument will not be used if the data have been prepared in advance via the `prep_data()` function.
#' @param prep_data Logical: Should the data be prepared before fitting the model? (Or have the data been prepared in advance by first running the `prep_data()` and `prep_data_cv()` functions)? If so, set `prep_data = FALSE`.) Defaults to `prep_data = TRUE`.
#' @param data A list of data produced by `prep_data()` followed by `prep_data_cv()`.
#' @param prefs An N × J matrix of numerical stimulus ratings or preference scores. These data are only required by the HBAM_R and HBAM_R_MINI models and will be ignored when fitting other models.
#' @param group_id Integer vector of length N identifying which group each respondent belongs to. The supplied vector should range from 1 to the total number of groups in the data, and all integers between these numbers should be represented in the supplied data. These data are only required by models with "MULTI" in their name and will be ignored when fitting other models.
#' @param K An integer above 2, specifying the number of folds to use in the analysis. Defaults to 10.
#' @param chains A positive integer specifying the number of Markov chains to use per fold. Defaults to 2.
#' @param warmup A positive integer specifying the number of warmup (aka burn-in) iterations per chain. It defaults to 1000. The number of warmup iterations should be smaller than `iter`.
#' @param iter A positive integer specifying the number of iterations for each chain (including warmup). It defaults to 3000 as running fewer chains for longer is a more efficient way to obtain a certain number of draws (and cross-validation can be computationally expensive).
#' @param seed An integer passed on to `set.seed` before creating the folds to increase reproducibility and comparability. Defaults to 1 and only applies to fold-creation when the argument `prep_data` is `TRUE`. The supplied `seed` argument is also used to generate seeds for the sampling algorithm.
#' @param control A named list of parameters to control the sampler's behavior. See the documentation for `rstan::stan` for more details.
#' @param sigma_alpha A positive numeric value specifying the standard deviation of the prior on the shift parameters in the FBAM model, or the standard deviation of the parameters' deviation from the group-means in FBAM_MULTI models. (This argument will be ignored by HBAM models.) Defaults to B / 4, where B measures the length of the survey scale as the number of possible placements on one side of the center.
#' @param sigma_beta A positive numeric value specifying the standard deviation of the prior on the logged stretch parameters in the FBAM model, or the standard deviation of the logged parameters' deviation from the group-means in FBAM_MULTI models. (This argument will be ignored by HBAM models.) Defaults to .35.
#' @param sigma_mu_alpha A positive numeric value specifying the standard deviation of the prior on the group-means of the shift parameters in MULTI-type models. Defaults to B / 5.
#' @param sigma_mu_beta A positive numeric value specifying the standard deviation of the prior on the group-means of the logged stretch parameters in MULTI-type models. Defaults to .3.
#' @param ... Arguments passed to `rstan::sampling()`.
#' @return A list of classes `kfold` and `loo`, which contains the following named elements:
#'    * `"estimates"`: A `1x2` matrix containing the ELPD estimate and its standard error. The columns have names `"Estimate"` and `"SE"`.
#'    * `"pointwise"`: A `Nx1` matrix with column name `"elpd_kfold"` containing the pointwise contributions for each data point.
#' @examples
#' \donttest{
#' # Loading and re-coding ANES 1980 data:
#' data(LC1980)
#' LC1980[LC1980 == 0 | LC1980 == 8 | LC1980 == 9] <- NA
#'
#' # Making a small subset of the data for illustration:
#' self <- LC1980[1:50, 1]
#' stimuli <- LC1980[1:50, -1]
#'
#' # Preparing to run chains in parallel using 2 cores via the future package:
#'   # Note: You would normally want to use all physical cores for this.
#' future::plan(future::multisession, workers = 2)
#'
#' # Performing 4-fold cross-validation for the HBAM_MINI model:
#'   # Note: You would typically want to run the chains for more iterations.
#' cv_hbam_mini <- hbam_cv(self, stimuli, model = "HBAM_MINI", K = 4,
#'                         chains = 1, warmup = 500, iter = 1000)
#'
#' # Performing 4-fold cross-validation for the FBAM model:
#' cv_fbam <- hbam_cv(self, stimuli, model = "FBAM", K = 4,
#'                         chains = 1, warmup = 500, iter = 1000)
#'
#' # Comparing the results using the loo package:
#' loo::loo_compare(list(HBAM_MINI = cv_hbam_mini,
#'                  FBAM = cv_fbam))
#'
#' # Stop the cluster of parallel sessions:
#' future::plan(future::sequential)
#' }

hbam_cv <- function(self = NULL, stimuli = NULL, model = "HBAM",
                    allow_miss = 0, req_valid = NA, req_unique = 2,
                    prefs = NULL, group_id = NULL, prep_data = TRUE, data = NULL, K = 10,
                    chains = 2,
                    warmup = 1000, iter = 3000,
                    seed = 1,
                    control = list(max_treedepth = 7),
                    sigma_alpha = NULL, sigma_beta = .35,
                    sigma_mu_alpha = NULL, sigma_mu_beta = .3, ...){

  logColMeansExp <- function(x) {
    S <- nrow(x)
    matrixStats::colLogSumExps(x) - log(S) # This is an alternative way of calculating the col means to taking sum of exps and dividing by N before then logging
  }

  if (prep_data == TRUE) {
    dat <- hbamr::prep_data(self = self, stimuli = stimuli, prefs = prefs, allow_miss = allow_miss, req_valid = req_valid, req_unique = req_unique, group_id = group_id)
    dat_l <- hbamr::prep_data_cv(data = dat, K = K, seed = seed)
  } else {
    dat_l <- data
  }

  n_fold <- length(dat_l)

  if (is.null(sigma_alpha)) { sigma_alpha <- dat_l[[1]]$B / 4.0 }
  if (is.null(sigma_mu_alpha)) { sigma_mu_alpha <- dat_l[[1]]$B / 5.0 }
  for (k in 1:n_fold) {
    dat_l[[k]]$sigma_alpha_fixed <- sigma_alpha
    dat_l[[k]]$sigma_mu_alpha <- sigma_mu_alpha
    dat_l[[k]]$sigma_beta_fixed <- sigma_beta
    dat_l[[k]]$sigma_mu_beta <- sigma_mu_beta
    dat_l[[k]] <- add_model_features(dat_l[[k]], model)
  }

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  # Function to run all chains and store only the log-likelihoods for held-out data to save memory:
  get_ll_list <- function() {
    p <- progressr::progressor(n_fold * chains)
    future.apply::future_lapply(1:(n_fold * chains),
                                function(i){
                                  k <- ceiling(i / chains) # Fold number
                                  init_l <- list(inits_omni(chain_id = i, dat = dat_l[[k]]))
                                  # Obtain chain:
                                  s <- rstan::sampling(stanmodels[[1]], data = dat_l[[k]], init = init_l, pars = "log_lik",
                                                       chains = 1, cores = 1, warmup = warmup, iter = iter, chain_id = i, seed = seed + i, control = control, refresh = 0, ...)
                                  # Calculate expected value of log-likelihood for each held-out observation:
                                  log_lik <- loo::extract_log_lik(s)
                                  draws <- dim(log_lik)[1]
                                  N_obs <- dim(log_lik)[2]
                                  heldout <- matrix(rep(dat_l[[k]]$holdout, each = draws), nrow = draws)
                                  log_lik_heldout <- log_lik * NA
                                  log_lik_heldout[heldout == 1] <- log_lik[heldout == 1]
                                  pointwise <-  matrix(logColMeansExp(log_lik_heldout), ncol= 1)

                                  p()

                                  return(pointwise)
                                },
                                future.seed = .Random.seed)

  }

  # Run all chains and show progress bar:
  pointwise_ll_list <- progressr::with_progress(get_ll_list())

  # Combine log-likelihoods from different chains by averaging, then combine values from different folds in a common vector:
  pointwise_ll <- vector()
  for(k in 1:n_fold){
    inchains <- (chains * k - (chains - 1)):(chains * k)
    pointwise_ll_mat <- t(matrix(unlist(pointwise_ll_list[inchains]), ncol = chains))
    pointwise_ll[dat_l[[k]]$holdout == 1] <- logColMeansExp(pointwise_ll_mat)[dat_l[[k]]$holdout == 1]
  }

  #out <- data.frame(
  #  ELPD = round(sum(pointwise_ll, na.rm = T), digits = 1),
  #  SE = round(as.numeric(sqrt(length(pointwise_ll) * var(pointwise_ll))), digits = 1))
  #rownames(out) <- model

  estimates = matrix(
    c(round(sum(pointwise_ll, na.rm = T), digits = 1),
      round(as.numeric(sqrt(length(pointwise_ll) * var(pointwise_ll))), digits = 1)),
    nrow = 1)
  rownames(estimates) <- "elpd_kfold"
  colnames(estimates) <- c("Estimate", "SE")

  pointwise_ll <- matrix(as.numeric(pointwise_ll), nrow = length(as.numeric(pointwise_ll)), ncol = 1)
  colnames(pointwise_ll) <- "elpd_kfold"

  out <- list("estimates" = estimates,
              "pointwise" = pointwise_ll)

  class(out) <- c("kfold", "loo")

  return(out)
}

