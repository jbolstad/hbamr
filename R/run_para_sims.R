run_para_sims <- function(nreps = 1, J_vals = 4, psi_vals = .89, tau_vals = .65,
                          it = 4000, warm = 1000, draws = 1500, ...) {

  # Create vectors containing settings for all simulations to be run
  vecs <- data.frame(
    J = rep(rep(J_vals, length(psi_vals) * length(tau_vals)), nreps),
    psi = rep(rep(rep(psi_vals, length(tau_vals)), each = length(J_vals)), nreps),
    tau = rep(rep(tau_vals, each = length(J_vals) * length(psi_vals)), nreps)
  )
  vecs$seed <- 1:nrow(vecs)

  # Run simulations in parallel
  res_list <-
    pbmcapply::pbmclapply(1:nrow(vecs), mc.cores = parallel::detectCores(logical = FALSE),
             function(i){
               results <- try(hbam:::run_sim(tau = vecs$tau[i], J = vecs$J[i], psi = vecs$psi[i],
                                  chains = 1, it = it, warm = warm, draws = draws,
                                  seed = vecs$seed[i], ...))
               if(inherits(results, "try-error")) {
                 error <- as.character(results)
                 results <- data.frame(rbind(rep(NA, 38)), vecs$seed[i], error)
                 names(results) <- c("N", "J", "nu", "tau", "psi", "cor_mean_spos", "ppos_mod1", "ppos_mod2", "ppos_mod3",
                   "cor_V_true_chi", "vpos_mod1", "vpos_mod2", "vpos_mod3","vpos_mean_mod1", "vpos_mean_mod2", "vpos_mean_mod3",
                   "mod3_median_psi", "mod3_mean_kappa","chi_cov_mod1", "chi_cov_mod1_ex_extr", "chi_cov_mod2", "chi_cov_mod3",
                   "extr_CI_mod1", "extr_CI_mod2", "extr_CI_mod3","rhat_mod1", "rhat_mod2", "rhat_mod3",
                   "rhat_all_mod1", "rhat_all_mod2", "rhat_all_mod3","neff_mod1", "neff_mod2", "neff_mod3","time_mod1", "time_mod2", "time_mod3", "time_it",
                   "seed", "error")
               }
               return(results)
             })

  results_full <- data.table::rbindlist(res_list)
  return(results_full)
}
