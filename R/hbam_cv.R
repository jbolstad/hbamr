#' Perform K-fold cross-validation of an HBAM-type model
#'
#' This function performs a K-fold cross-validation of an HBAM-type model in order to estimate the expected log pointwise predictive density for a new dataset (ELPD).
#'
#' @export
#' @param self A numerical vector of N ideological self-placements. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param stimuli An N × J matrix of numerical stimulus placements, where J is the number of stimuli. Any missing data must be coded as NA. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param model Character: Name of the model to be used. Defaults to `"HBAM"`.
#' @param allow_miss Integer specifying how many missing stimulus positions to be accepted for an individual still to be included in the analysis. This argument will not be used if the data have been prepared in advance via the `prep_data` function. Defaults to 0.
#' @param req_valid Integer specifying how many valid observations to require for a respondent to be included in the analysis. The default is `req_valid = J - allow_miss`, but if specified, `req_valid` takes precedence. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param req_unique Integer specifying how may unique positions on the ideological scale each respondent is required to have used when placing the stimuli in order to be included in the analysis. The default is `req_unique = 2`. This argument will not be used if the data have been prepared in advance via the `prep_data` function.
#' @param prep_data Logical: Should the data be prepared before fitting the model? (Or have the data been prepared in advance by first running the `prep_data` and `prep_data_cv` functions)? If so, set `prep_data = FALSE`.) Defaults to `prep_data = TRUE`.
#' @param data A list of data produced by `prep_data` followed by `prep_data_cv`.
#' @param prefs An N × J matrix of numerical stimulus ratings or preference scores. These data are only required by the `"HBAM_R"` and `"HBAM_R_MINI"` models and will be ignored when fitting other models.
#' @param K An integer above 2, specifying the number of folds to use in the analysis. Defaults to 10.
#' @param chains A positive integer specifying the number of Markov chains to use for each model fit. Defaults to 2.
#' @param cores The number of cores to use when executing the Markov chains in parallel. Defaults to `parallel::detectCores(logical = FALSE)`. The function is parallelized so that users can specify a higher number of cores than chains and run chains for different folds simultaneously to save time.
#' @param warmup A positive integer specifying the number of warmup (aka burn-in) iterations per chain. If step-size adaptation is on (which it is by default), this also controls the number of iterations for which adaptation is run (and hence these warmup samples should not be used for inference). The number of warmup iterations should be smaller than `iter`.
#' @param iter A positive integer specifying the number of iterations for each chain (including warmup).
#' @param thin A positive integer specifying the period for saving samples.
#' @param control A named list of parameters to control the sampler's behavior. See the details in the documentation for the control argument in the `stan` function in the `rstan` package.
#' @param seed An integer passed on to `set.seed` before creating the folds to increase reproducibility and comparability. Defaults to 1 and only applies to fold-creation when the argument `prep_data` is `TRUE`. The supplied `seed` argument is also used to generate seeds for the sampling algorithm.
#' @return A data frame containing the estimated ELPD and its standard error.
#' @examples
#' \dontrun{
#' # Loading and re-coding ANES 1980 data:
#' data(LC1980)
#' dat <- LC1980
#' dat[dat == 0 | dat == 8 | dat == 9] <- NA
#' self <- dat[, 1]
#' stimuli <- dat[, -1]
#'
#' # Performing 10-fold cross-validation for the HBAM model:
#' cv_hbam <- hbam_cv(self, stimuli, model = "HBAM")
#' }

hbam_cv <- function(self = NULL, stimuli = NULL, model = "HBAM",
                    allow_miss = 0, req_valid = NA, req_unique = 2,
                    prefs = NULL, prep_data = TRUE, data = NULL, K = 10,
                    chains = 2, cores = parallel::detectCores(logical = FALSE),
                    warmup = 1000, iter = 4000,
                    thin = 1,
                    control = list(adapt_delta = .6), seed = 1){

  logColMeansExp <- function(x) {
    S <- nrow(x)
    matrixStats::colLogSumExps(x) - log(S) # This is an alternative way of calculating the col means to taking sum of exps and dividing by N before then logging
  }

  if (prep_data == TRUE) {
    dat <- hbamr::prep_data(self = self, stimuli = stimuli, prefs = prefs, allow_miss = allow_miss, req_valid = req_valid, req_unique = req_unique)
    dat_l <- hbamr::prep_data_cv(data = dat, K = K, seed = seed)
  } else {
    dat_l <- data
  }

  n_fold <- length(dat_l)
  # Run all chains in parallel and store only the log-likelihoods for held-out data to save memory:
  pointwise_ll_list <-
    pbmcapply::pbmclapply(1:(n_fold * chains),
                          function(i){
                            k <- ceiling(i / chains) # Fold number
                            # Generate inits
                            set.seed(seed + i)
                            init_l <- list(inits[[model]](chain_id = i, dat = dat_l[[k]]))
                            # Obtain chain:
                            s <- rstan::sampling(stanmodels[[model]], data = dat_l[[k]], init = init_l,
                                                 chains = 1, cores = 1, warmup = warmup, iter = iter, thin = thin, control = control, chain_id = i, seed = seed + i)
                            # Calculate expected value of log-likelihood for each held-out observation:
                            log_lik <- loo::extract_log_lik(s)
                            draws <- dim(log_lik)[1]
                            N_obs <- dim(log_lik)[2]
                            heldout <- matrix(rep(dat_l[[k]]$holdout, each = draws), nrow = draws)
                            log_lik_heldout <- log_lik * NA
                            log_lik_heldout[heldout == 1] <- log_lik[heldout == 1]
                            pointwise <-  matrix(logColMeansExp(log_lik_heldout), ncol= 1)

                            return(pointwise)
                          },
                          mc.cores = cores)

  # Combine log-likelihoods from different chains by averaging, then combine values from different folds in a common vector:
  pointwise_ll <- vector()
  for(k in 1:n_fold){
    inchains <- (chains * k - (chains - 1)):(chains * k)
    pointwise_ll_mat <- t(matrix(unlist(pointwise_ll_list[inchains]), ncol = chains))
    pointwise_ll[dat_l[[k]]$holdout == 1] <- logColMeansExp(pointwise_ll_mat)[dat_l[[k]]$holdout == 1]
  }

  out <- data.frame(
    ELPD = round(sum(pointwise_ll, na.rm = T), digits = 1),
    SE = round(as.numeric(sqrt(length(pointwise_ll) * var(pointwise_ll))), digits = 1))
  rownames(out) <- model

  return(out)
}
