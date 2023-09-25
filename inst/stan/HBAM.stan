data {
  int<lower = 1> N;                       // n of individuals
  int<lower = 1> J;                       // n of items
  int<lower = 1> N_obs;                   // n of observations
  array[N_obs] int<lower = 1> ii;         // index i in matrix
  array[N_obs] int<lower = 1> jj;         // index j in matrix
  int<lower = 1> B;                       // length of scale -1 / 2
  int<lower = 1, upper = J> L;            // left pole
  int<lower = 1, upper = J> R;            // right pole
  array[N_obs] int<lower = -B, upper = B> Y; // reported stimuli positions
  vector<lower = -B, upper = B>[N] V;     // reported self-placements
  int<lower=0, upper=1> CV;               // indicator of cross-validation
  array[N_obs] int<lower=0, upper=1> holdout; // holdout for cross-validation
}

transformed data {
  real<lower = 0> sigma_chi_prior_rate = (5 - 1) / (B / 2.0);
  real<lower = 0> sigma_alpha_prior_rate = (2 - 1) / (B / 5.0);
  real<lower = 0> tau_prior_rate = (2 - 1) / (B / 5.0);
}

parameters {
  matrix[N, 2] alpha_raw;                 // shift parameter, split, raw
  matrix[N, 2] beta_raw;                  // stretch parameter, split, raw
  ordered[2] theta_lr;                    // left and right pole
  array[J] real theta_raw;                // remaining stimuli
  real<lower = 0> sigma_alpha;            // sd of alpha
  real<lower = 0, upper = 2> sigma_beta;  // sd of log(beta)
  real<lower = 0> sigma_chi;              // sd of chi0
  real<lower = 3, upper = 30> nu;         // concentration of etas
  real<lower = 0> tau;                    // scale of errors
  vector<lower = 0>[N] eta;               // mean ind. error variance x J^2
  simplex[J] rho;                         // stimuli-shares of variance
  vector<lower = 0, upper = 1>[N] lambda; // mixing proportion, flipping
  real<lower = .5, upper = 1> psi;        // mean of prior on lambda
  real<lower = 2, upper = 100> delta;     // concentration of prior on lambda
  matrix[N, 2] chi0;                      // latent respondent positions, split
}

transformed parameters {
  real<lower=0> alpha_lambda = delta * psi; // reparameterization
  real<lower=0> beta_lambda = delta * (1 - psi);
  array[J] real theta;                    // latent stimuli position
  matrix[N, 2] alpha0;                    // shift parameter, split
  matrix[N, 2] beta0;                     // stretch parameter, split
  vector[N_obs] log_lik;                  // pointwise log-likelihood for Y
  vector[N] log_lik_V;                    // pointwise log-likelihood for V
  real<lower = 0> eta_scale = tau * J;
  real<lower=0> min_rho = min(rho);
  theta = theta_raw;
  theta[L] = theta_lr[1];                 // safeguard to ensure identification
  theta[R] = theta_lr[2];
  alpha0[, 1] = alpha_raw[, 1] * sigma_alpha; // non-centered specifications
  alpha0[, 2] = alpha_raw[, 2] * sigma_alpha;
  beta0[, 1] = exp(beta_raw[, 1] * sigma_beta);
  beta0[, 2] = -exp(beta_raw[, 2] * sigma_beta);

  for (n in 1:N_obs) {
    log_lik[n] = log_mix( lambda[ii[n]],
      normal_lpdf(Y[n] | alpha0[ii[n], 1] + beta0[ii[n], 1] * theta[jj[n]],
        sqrt(eta[ii[n]]) * rho[jj[n]]),
      normal_lpdf(Y[n] | alpha0[ii[n], 2] + beta0[ii[n], 2] * theta[jj[n]],
        sqrt(eta[ii[n]]) * rho[jj[n]]) );
  }
  for (i in 1:N) {
    log_lik_V[i] = log_mix( lambda[i],
      normal_lpdf(V[i] | alpha0[i, 1] + beta0[i, 1] * chi0[i, 1],
        sqrt(eta[i]) * min_rho),
      normal_lpdf(V[i] | alpha0[i, 2] + beta0[i, 2] * chi0[i, 2],
        sqrt(eta[i]) * min_rho) );
  }
}

model {
  theta_raw ~ normal(0, B);
  theta_lr ~ normal(0, B);
  alpha_raw[, 1] ~ normal(0, 1);
  alpha_raw[, 2] ~ normal(0, 1);
  sigma_alpha ~ gamma(2, sigma_alpha_prior_rate);
  beta_raw[, 1] ~ normal(0, 1);
  beta_raw[, 2] ~ normal(0, 1);
  sigma_beta ~ gamma(3, 10);
  chi0[, 1] ~ normal(0, sigma_chi);
  chi0[, 2] ~ normal(0, sigma_chi);
  sigma_chi ~ gamma(5, sigma_chi_prior_rate);
  eta ~ scaled_inv_chi_square(nu, eta_scale);
  nu ~ gamma(25, 2.5);
  tau ~ gamma(2, tau_prior_rate);
  rho ~ dirichlet(rep_vector(5, J));
  lambda ~ beta(alpha_lambda, beta_lambda);
  psi ~ beta(8.5, 1.5);
  delta - 2 ~ gamma(2, .1);

  target += sum(log_lik_V);
  if(CV == 0)
    target += sum(log_lik);
  else
    for (n in 1:N_obs) {
      if(holdout[n] == 0)
        target += log_lik[n];
    }
}

generated quantities {
  vector[N] kappa = to_vector(bernoulli_rng(lambda));
  vector[N] chi = (kappa .* chi0[, 1]) + ((1 - kappa) .* chi0[, 2]);
  vector[N] alpha = (kappa .* alpha0[, 1]) + ((1 - kappa) .* alpha0[, 2]);
  vector[N] beta = (kappa .* beta0[, 1]) + ((1 - kappa) .* beta0[, 2]);
}
