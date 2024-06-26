data {
  int<lower = 1> N;                       // n of individuals
  int<lower = 1> J;                       // n of items
  int<lower = 1> G;                       // n of groups
  int<lower = 1> N_obs;                   // n of observations
  array[N_obs] int<lower = 1> ii;         // index i in matrix
  array[N_obs] int<lower = 1> jj;         // index j in matrix
  array[N] int<lower = 1> gg;             // group-index of i
  int<lower = 1> B;                       // length of scale -1 / 2
  int<lower = 1, upper = J> L;            // left pole
  int<lower = 1, upper = J> R;            // right pole
  vector<lower = -B, upper = B>[N_obs] Y; // reported stimuli positions
  vector<lower = -B, upper = B>[N] V;     // reported self-placements
  int<lower = 0, upper = 1> CV;           // indicator of cross-validation
  vector<lower = 0, upper = 1>[N_obs] holdout; // holdout for cross-validation
  real<lower = 0> sigma_alpha;            // sd of prior on alpha
  real<lower = 0> sigma_beta;             // sd of prior on log(beta)
  real<lower = 0> sigma_mu_alpha;         // sd of prior on mu_alpha
  real<lower = 0> sigma_mu_beta;          // sd of prior on mu_beta
  int<lower = 0, upper = 1> MCMC;         // indicator of fitting method
}

transformed data {
  real nu = 6;                            // concentration of etas
  real tau = B / 4.0;                     // scale of prior on errors
  real eta_scale = tau * J;
  vector<lower = 0, upper = 1>[N_obs] not_holdout = 1 - holdout;
  real mean_mu_simplexes = 1.0 / G;       // for later scaling of simplexes
  real sd_mu_simplexes = sqrt(mean_mu_simplexes * (1 - mean_mu_simplexes) / (50 * G + 1));
}

parameters {
  vector[N] alpha_raw;                    // shift parameter, raw
  vector[N] beta_raw;                     // stretch parameter, raw
  ordered[2] theta_lr;                    // left and right pole
  array[J] real theta_raw;                // remaining stimuli
  simplex[G] mu_alpha_raw;                // group-level mean of alpha, raw
  simplex[G] mu_beta_raw;                 // group-level mean of log(beta), raw
  vector<lower = 0>[N] eta;               // mean ind. error variance x J^2
  simplex[J] rho;                         // stimuli-shares of variance
}

transformed parameters {
  array[J] real theta;                    // latent stimuli position
  vector[N] alpha;                        // shift parameter
  vector[N] beta;                         // stretch parameter
  vector[N_obs] log_lik;                  // pointwise log-likelihood for Y
  vector[G] mu_alpha = ((mu_alpha_raw - mean_mu_simplexes) / sd_mu_simplexes) * sigma_mu_alpha;
  vector[G] mu_beta = ((mu_beta_raw - mean_mu_simplexes) / sd_mu_simplexes) * sigma_mu_beta;
  theta = theta_raw;
  theta[L] = theta_lr[1];                 // safeguard to ensure identification
  theta[R] = theta_lr[2];

  for (i in 1:N) {
    alpha[i] = alpha_raw[i] * sigma_alpha + mu_alpha[gg[i]]; // non-centered specifications
    beta[i] = exp(beta_raw[i] * sigma_beta + mu_beta[gg[i]]);
  }

  for (n in 1:N_obs) {
    log_lik[n] = normal_lpdf(Y[n] | alpha[ii[n]] + beta[ii[n]] * theta[jj[n]],
        sqrt(eta[ii[n]]) * rho[jj[n]]);
  }

}

model {
  theta_raw ~ normal(0, B / 2.0);
  theta_lr ~ normal(0, B / 2.0);
  alpha_raw ~ normal(0, 1);
  beta_raw ~ normal(0, 1);
  mu_alpha_raw ~ dirichlet(rep_vector(50, G));
  mu_beta_raw ~ dirichlet(rep_vector(50, G));
  eta ~ scaled_inv_chi_square(nu, eta_scale);
  rho ~ dirichlet(rep_vector(50, J));

  if (CV == 0)
    target += sum(log_lik);
  else
    target += sum(log_lik .* not_holdout);
}

generated quantities {
  vector[N] chi;
  real<lower = 0> min_rho = min(rho);
  if (MCMC == 1)
    chi = (V - to_vector(normal_rng(0, sqrt(eta) * min_rho)) - alpha) ./ beta;
  else
    chi = (V - alpha) ./ beta;
}
