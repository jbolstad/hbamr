# Function to run simulations:
run_sim <- function(N = 500, J = 4, nu = 8.4, tau = .57, psi = .89,
                    chains = 4, it = 2000, warm = 1000, draws = 1000, seed = 1, ...) {

  ptm.it <- proc.time()
  dat <- hbamr:::create_data(N = N, J = J, nu = nu, tau = tau, psi = psi, seed = seed, ...)
  thin <- round((it - warm) / (draws / chains))

  ptm <- proc.time()
  mod1_fit <- hbamr::hbam(data = dat, model = "BAM", prep_data = FALSE,
                   chains = chains, cores = chains, refresh = 0,
                   warmup = warm, iter = it, thin = thin, seed = seed)
  t_mod1 <- (proc.time() - ptm)[[3]]

  ptm <- proc.time()
  mod2_fit <- hbamr::hbam(data = dat, model = "HBAM_0", prep_data = FALSE,
                   chains = chains, cores = chains, refresh = 0,
                   warmup = warm, iter = it, thin = thin, seed = seed)
  t_mod2 <- (proc.time() - ptm)[[3]]

  ptm <- proc.time()
  mod3_fit <- hbamr::hbam(data = dat, model = "HBAM", prep_data = FALSE,
                   chains = chains, cores = chains, refresh = 0,
                   warmup = warm, iter = it, thin = thin, seed = seed)
  t_mod3 <- (proc.time() - ptm)[[3]]

  mod1_theta <- rstan::summary(mod1_fit, pars = "theta")[[1]][, 6]
  mod1_chi <- rstan::summary(mod1_fit, pars = "chi")[[1]][, 6]
  mod1_mean_chi <- rstan::summary(mod1_fit, pars = "chi")[[1]][, 1]

  mod2_theta <- rstan::summary(mod2_fit, pars = "theta")[[1]][, 6]
  mod2_chi <- rstan::summary(mod2_fit, pars = "chi")[[1]][, 6]
  mod2_mean_chi <- rstan::summary(mod2_fit, pars = "chi")[[1]][, 1]

  mod3_theta <- rstan::summary(mod3_fit, pars = "theta")[[1]][, 6]
  mod3_chi <- rstan::summary(mod3_fit, pars = "chi")[[1]][, 6]
  mod3_mean_chi <- rstan::summary(mod3_fit, pars = "chi")[[1]][, 1]

  mod3_median_psi <- rstan::summary(mod3_fit, pars = "psi")[[1]][6]
  mod3_mean_kappa <- mean(rstan::summary(mod3_fit, pars = "kappa")[[1]][, 1])

  mod1_extr <- hbamr:::extreme_CIs(mod1_fit)
  mod1_cov <- hbamr:::CI_coverage(mod1_fit, dat)

  mod2_extr <- hbamr:::extreme_CIs(mod2_fit)
  mod2_cov <- hbamr:::CI_coverage(mod2_fit, dat)

  mod3_extr <- hbamr:::extreme_CIs(mod3_fit)
  mod3_cov <- hbamr:::CI_coverage(mod3_fit, dat)

  results <- data.frame(N, J, nu, tau, psi, cor(dat$mean_spos, dat$true_theta),
    cor(mod1_theta, dat$true_theta), cor(mod2_theta, dat$true_theta), cor(mod3_theta, dat$true_theta),
    cor(dat$V, dat$true_chi), cor(mod1_chi, dat$true_chi), cor(mod2_chi, dat$true_chi), cor(mod3_chi, dat$true_chi),
    cor(mod1_mean_chi, dat$true_chi), cor(mod2_mean_chi, dat$true_chi), cor(mod3_mean_chi, dat$true_chi),
    mod3_median_psi, mod3_mean_kappa,
    mean(mod1_cov), mean(mod1_cov[!mod1_extr]), mean(mod2_cov), mean(mod3_cov),
    mean(mod1_extr), mean(mod2_extr), mean(mod3_extr),
    max(rstan::summary(mod1_fit, pars = c("theta", "chi"))$summary[, 10]),
    max(rstan::summary(mod2_fit, pars = c("theta", "chi"))$summary[, 10]),
    max(rstan::summary(mod3_fit, pars = c("theta", "chi"))$summary[, 10]),
    max(rstan::summary(mod1_fit)$summary[, 10]),
    max(rstan::summary(mod2_fit)$summary[, 10]),
    max(rstan::summary(mod3_fit)$summary[, 10]),
    min(rstan::summary(mod1_fit, pars = c("theta", "chi"))$summary[, 9]),
    min(rstan::summary(mod2_fit, pars = c("theta", "chi"))$summary[, 9]),
    min(rstan::summary(mod3_fit, pars = c("theta", "chi"))$summary[, 9]),
    t_mod1, t_mod2, t_mod3, (proc.time() - ptm.it)[[3]], seed, "no")

  names(results) <- c(
    "N", "J", "nu", "tau", "psi", "cor_mean_spos",
    "ppos_mod1", "ppos_mod2", "ppos_mod3",
    "cor_V_true_chi", "vpos_mod1", "vpos_mod2", "vpos_mod3",
    "vpos_mean_mod1", "vpos_mean_mod2", "vpos_mean_mod3",
    "mod3_median_psi", "mod3_mean_kappa",
    "chi_cov_mod1", "chi_cov_mod1_ex_extr", "chi_cov_mod2", "chi_cov_mod3",
    "extr_CI_mod1", "extr_CI_mod2", "extr_CI_mod3",
    "rhat_mod1", "rhat_mod2", "rhat_mod3",
    "rhat_all_mod1", "rhat_all_mod2", "rhat_all_mod3",
    "neff_mod1", "neff_mod2", "neff_mod3",
    "time_mod1", "time_mod2", "time_mod3", "time_it",
    "seed", "error")

  return(results)
}
