data {
  int<lower = 1> N;                       // n of individuals
  int<lower = 1> J;                       // n of items
  int<lower = 1> N_obs;                   // n of observations
  array[N_obs] int<lower = 1> ii;         // index i in matrix
  array[N_obs] int<lower = 1> jj;         // index j in matrix
  int<lower = 1> B;                       // length of scale -1 / 2
  int<lower = 1, upper = J> L;            // left pole
  int<lower = 1, upper = J> R;            // right pole
  vector<lower = -B, upper = B>[N_obs] Y; // reported stimuli positions
  vector<lower = -B, upper = B>[N] V;     // reported self-placements
  int<lower = 0, upper = 1> CV;           // indicator of cross-validation
  vector<lower = 0, upper = 1>[N_obs] holdout; // holdout for cross-validation
  real<lower = 0> sigma_alpha;            // sd of prior on alpha
  real<lower = 0> sigma_beta;             // sd of prior on log(beta)
  int<lower = 0, upper = 1> MCMC;         // indicator of fitting method
}

transformed data {
  real<lower = 0> tau_prior_rate = (2 - 1) / (B / 5.0);
  vector<lower = 0, upper = 1>[N_obs] not_holdout = 1 - holdout;
}

parameters {
  matrix[N, 2] alpha_raw;                 // shift parameter, split, raw
  matrix[N, 2] beta_raw;                  // stretch parameter, split, raw
  ordered[2] theta_lr;                    // left and right pole
  array[J] real theta_raw;                // remaining stimuli
  real<lower = 0> tau;                    // scale of errors
  vector<lower = 0, upper = 1>[N] lambda; // mixing proportion, flipping
}

transformed parameters {
  array[J] real theta;                    // latent stimuli position
  matrix[N, 2] alpha0;                    // shift parameter, split
  matrix[N, 2] beta0;                     // stretch parameter, split
  vector[N_obs] log_lik;                  // pointwise log-likelihood for Y
  theta = theta_raw;
  theta[L] = theta_lr[1];                 // safeguard to ensure identification
  theta[R] = theta_lr[2];
  alpha0[, 1] = alpha_raw[, 1] * sigma_alpha; // non-centered specifications
  alpha0[, 2] = alpha_raw[, 2] * sigma_alpha;
  beta0[, 1] = exp(beta_raw[, 1] * sigma_beta);
  beta0[, 2] = -exp(beta_raw[, 2] * sigma_beta);

  for (n in 1:N_obs) {
    log_lik[n] = log_mix( lambda[ii[n]],
      normal_lpdf(Y[n] | alpha0[ii[n], 1] + beta0[ii[n], 1] * theta[jj[n]], tau),
      normal_lpdf(Y[n] | alpha0[ii[n], 2] + beta0[ii[n], 2] * theta[jj[n]], tau) );
  }
}

model {
  theta_raw ~ normal(0, B);
  theta_lr ~ normal(0, B);
  alpha_raw[, 1] ~ normal(0, 1);
  alpha_raw[, 2] ~ normal(0, 1);
  beta_raw[, 1] ~ normal(0, 1);
  beta_raw[, 2] ~ normal(0, 1);
  tau ~ gamma(2, tau_prior_rate);
  lambda ~ beta(2, 1);

  if (CV == 0)
    target += sum(log_lik);
  else
    target += sum(log_lik .* not_holdout);
}

generated quantities {
  vector[N] kappa;
  vector[N] chi;
  if (MCMC == 1)
    kappa = to_vector(bernoulli_rng(lambda));
  else
    kappa = to_vector(round(lambda)); // Rounding to MAP instead of sampling
  vector[N] alpha = (kappa .* alpha0[, 1]) + ((1 - kappa) .* alpha0[, 2]);
  vector[N] beta = (kappa .* beta0[, 1]) + ((1 - kappa) .* beta0[, 2]);
  if (MCMC == 1)
    chi = (V - to_vector(normal_rng(0, rep_vector(tau, N))) - alpha) ./ beta;
  else
    chi = (V - alpha) ./ beta;
}
