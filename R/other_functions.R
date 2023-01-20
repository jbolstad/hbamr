get_post_mode_x <- function(d) {dd <- density(d); dd$x[which.max(dd$y)]}
get_post_mode_y <- function(d) {dd <- density(d); dd$y[which.max(dd$y)]}

extreme_CIs <- function(object, prob = c(.025, .975), par = "chi", check_fac = 1.5) {
  draws <- as.data.frame(rstan::extract(object, pars = "chi")$chi)
  draws <- tidyr::pivot_longer(draws, starts_with("V"), values_to = "chi")
  qtile <- quantile(draws$chi, prob)
  prange <- qtile[2] - qtile[1]
  CI <- rstan::summary(object, pars = par)[[1]][, c(4, 8)]
  ranges <- CI[, 2] - CI[, 1]
  rm(draws)
  return(ranges > prange * check_fac)
}

out_of_pop_dist <- function(object, prob = c(.025, .975), par = "chi", check_fac = 2) {
  mean_median<- rstan::summary(object, pars = par)[[1]][, c(1, 6)]
  draws <- as.data.frame(rstan::extract(object, pars = "chi")$chi)
  draws <- tidyr::pivot_longer(draws, starts_with("V"), values_to = "chi")
  median_all <- median(draws$chi)
  mad_all <- mad(draws$chi)
  range <- c(median_all - mad_all * check_fac, median_all + mad_all * check_fac)
  rm(draws)
  return(data.frame(mean_outof_range = !(mean_median[, 1] > range[1] & mean_median[, 1] < range[2]),
                    median_outof_range = !(mean_median[, 2] > range[1] & mean_median[, 2] < range[2])))
}

CI_coverage <- function(object, dat, par = "chi") {
  CI <- rstan::summary(object, pars = par)[[1]][, c(4, 8)]
  median_est <- rstan::summary(object, pars = par)[[1]][, 6]
  CI <- (CI - mean(median_est)) / sd(median_est)
  truevals <- dat[[paste("true_", par, sep = "")]]
  truevals <- (truevals - mean(truevals)) / sd(truevals)
  in_CI <- truevals > CI[, 1] & truevals < CI[, 2]
  return(in_CI)
}

