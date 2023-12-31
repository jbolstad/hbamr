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
  array[N_obs] int<lower = -B, upper = B> Y; // reported stimuli positions
  vector<lower = -B, upper = B>[N] V;     // reported self-placements
  int<lower=0, upper=1> CV;               // indicator of cross-validation
  array[N_obs] int<lower=0, upper=1> holdout; // holdout for cross-validation
}

transformed data {
  real<lower = 0> sigma_alpha_prior_rate = (2 - 1) / (B / 5.0);
  real<lower = 0> tau_prior_rate = (2 - 1) / (B / 5.0);
  real<lower = 0> sigma_mu_alpha = B / 5.0; // sd of mu_alpha
  real sigma_mu_beta = .25;               // sd of mu_beta
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
  real<lower = 0> sigma_alpha;            // sd of alpha
  real<lower = 0, upper = 2> sigma_beta;  // sd of log(beta)
  real<lower = 3, upper = 30> nu;         // concentration of etas
  real<lower = 0> tau;                    // scale of errors
  vector<lower = 0>[N] eta;               // mean ind. error variance x J^2
  simplex[J] rho;                         // stimuli-shares of variance
}

transformed parameters {
  vector[N] alpha;                        // shift parameter
  vector[N] beta;                         // stretch parameter
  vector[N] chi;                          // latent respondent positions
  array[J] real theta;                    // latent stimuli position
  vector[N_obs] log_lik;                  // pointwise log-likelihood for Y
  vector[G] mu_alpha = ((mu_alpha_raw - mean_mu_simplexes) / sd_mu_simplexes) * sigma_mu_alpha;
  vector[G] mu_beta = ((mu_beta_raw - mean_mu_simplexes) / sd_mu_simplexes) * sigma_mu_beta;
  real<lower = 0> eta_scale = tau * J;
  theta = theta_raw;
  theta[L] = theta_lr[1];                 // safeguard to ensure identification
  theta[R] = theta_lr[2];

  for (i in 1:N) {
    alpha[i] = alpha_raw[i] * sigma_alpha + mu_alpha[gg[i]]; // non-centered specifications
    beta[i] = exp(beta_raw[i] * sigma_beta + mu_beta[gg[i]]);
  }
  chi = ((V - alpha) ./ beta);

  for (n in 1:N_obs) {
    log_lik[n] = normal_lpdf(Y[n] | alpha[ii[n]] + beta[ii[n]] * theta[jj[n]],
      sqrt(eta[ii[n]]) * rho[jj[n]]);
  }
}

model {
  theta_raw ~ normal(0, B);
  theta_lr ~ normal(0, B);
  alpha_raw ~ normal(0, 1);
  sigma_alpha ~ gamma(2, sigma_alpha_prior_rate);
  beta_raw ~ normal(0, 1);
  sigma_beta ~ gamma(3, 10);
  mu_alpha_raw ~ dirichlet(rep_vector(50, G));
  mu_beta_raw ~ dirichlet(rep_vector(50, G));
  eta ~ scaled_inv_chi_square(nu, eta_scale);
  nu ~ gamma(25, 2.5);
  tau ~ gamma(2, tau_prior_rate);
  rho ~ dirichlet(rep_vector(5, J));

  if(CV == 0)
    target += sum(log_lik);
  else
    for (n in 1:N_obs) {
      if(holdout[n] == 0)
        target += log_lik[n];
    }
}
