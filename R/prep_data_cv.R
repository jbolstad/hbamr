#' Prepare data for a K-fold cross-validation of an HBAM model
#'
#' This function turns data prepared for `hbam` into a list of K versions, where each version includes a different vector identifying holdout-data.
#'
#' @export
#' @param data A list of data produced by `prep_data`.
#' @param K An integer above 2, specifying the number of folds to use in the analysis. Defaults to 10.
#' @param seed An integer passed on to `set.seed` before creating the folds to increase reproducibility. Defaults to 1.
#' @return A list of K data objects where each version includes a different vector identifying holdout-data.
#' @examples
#' # Loading and re-coding ANES 1980 data:
#' data(LC1980)
#' LC1980[LC1980 == 0 | LC1980 == 8 | LC1980 == 9] <- NA
#' self <- LC1980[, 1]
#' stimuli <- LC1980[, -1]
#' dat <- prep_data(self, stimuli)
#'
#' # Prepare data for cross-validation:
#' dat_cv <- prep_data_cv(dat, K = 10)
#'

prep_data_cv <- function(data, K = 10, seed = 1) {
  set.seed(seed)
  hh <- loo::kfold_split_stratified(K = K, x = data$ii)
  holdout_k <- matrix(0, nrow = data$N_obs, ncol = K)
  for(i in 1:data$N_obs) holdout_k[i, hh[i]] <- 1
  holdout_k <- split(holdout_k, rep(1:ncol(holdout_k), each = nrow(holdout_k)))
  data_l <- rep(list(data), K)
  for(i in 1:K) {
    data_l[[i]]$holdout <- holdout_k[[i]]
    data_l[[i]]$CV <- 1
  }
  return(data_l)
}
