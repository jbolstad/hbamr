data {
  int<lower = 1> N;                       // n of individuals
  int<lower = 1> J;                       // n of items
  int<lower = 1> N_obs;                   // n of observations
  int<lower = 1> ii[N_obs];               // index i in matrix
  int<lower = 1> jj[N_obs];               // index j in matrix
  int<lower = 1> B;                       // length of scale -1 / 2
  int<lower = -B, upper = B> Y[N_obs];    // reported stimuli positions
  vector<lower = -B, upper = B>[N] V;     // reported self-placements
  int<lower=0, upper=1> CV;               // indicator of cross-validation
  int<lower=0, upper=1> holdout[N_obs];   // holdout for cross-validation
}

transformed data {
  real<lower = 0> sigma_chi_prior_rate = (5 - 1) / (B / 2.0);
  real<lower = 0> sigma_alpha_prior_rate = (2 - 1) / (B / 5.0);
  real<lower = 0> tau_prior_rate = (2 - 1) / (B / 5.0);
}

parameters {
  vector[N] alpha_raw;                    // shift parameter, raw
  vector[N] beta_raw;                     // stretch parameter, raw
  real theta[J];                          // latent stimuli position
  real<lower = 0> sigma_alpha;            // sd of alpha
  real<lower = 0, upper = 2> sigma_beta;  // sd of log(beta)
  real<lower = 0> sigma_chi;              // sd of chi0
  real<lower = 3, upper = 30> nu;         // concentration of etas
  real<lower = 0> tau;                    // scale of errors
  vector<lower = 0>[N] eta;               // mean ind. error variance x J^2
  simplex[J] rho;                         // stimuli-shares of variance
  vector[N] chi;                          // latent respondent positions
}

transformed parameters {
  vector[N] alpha;                        // shift parameter
  vector[N] beta;                         // stretch parameter
  vector[N_obs] log_lik;                  // pointwise log-likelihood for Y
  vector[N] log_lik_V;                    // pointwise log-likelihood for V
  real<lower = 0> eta_scale = tau * J;
  real<lower=0> min_rho = min(rho);
  alpha = alpha_raw * sigma_alpha;        // non-centered specifications
  beta = exp(beta_raw * sigma_beta);

  for (n in 1:N_obs) {
    log_lik[n] = normal_lpdf(Y[n] | alpha[ii[n]] + beta[ii[n]] * theta[jj[n]],
      sqrt(eta[ii[n]]) * rho[jj[n]]);
  }
  for (i in 1:N) {
    log_lik_V[i] = normal_lpdf(V[i] | alpha[i] + beta[i] * chi[i], sqrt(eta[i]) * min_rho);
  }
}

model {
  theta ~ normal(0, B * 2);
  alpha_raw ~ normal(0, 1);
  sigma_alpha ~ gamma(2, sigma_alpha_prior_rate);
  beta_raw ~ normal(0, 1);
  sigma_beta ~ gamma(3, 10);
  chi ~ normal(0, sigma_chi);
  sigma_chi ~ gamma(5, sigma_chi_prior_rate);
  eta ~ scaled_inv_chi_square(nu, eta_scale);
  nu ~ gamma(25, 2.5);
  tau ~ gamma(2, tau_prior_rate);
  rho ~ dirichlet(rep_vector(5, J));

  target += sum(log_lik_V);
  if(CV == 0)
    target += sum(log_lik);
  else
    for (n in 1:N_obs) {
      if(holdout[n] == 0)
        target += log_lik[n];
    }
}
