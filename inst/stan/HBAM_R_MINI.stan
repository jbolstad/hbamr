data {
  int<lower = 1> N;                       // n of individuals
  int<lower = 1> J;                       // n of items
  int<lower = 1> N_obs;                   // n of observations
  array[N_obs] int<lower = 1> ii;               // index i in matrix
  array[N_obs] int<lower = 1> jj;               // index j in matrix
  int<lower = 1> B;                       // length of scale -1 / 2
  int<lower = 1, upper = J> L;            // left pole
  int<lower = 1, upper = J> R;            // right pole
  array[N_obs] int<lower = -B, upper = B> Y;    // reported stimuli positions
  array[N_obs] real<lower = 0, upper = 1> U;    // reported voter preferences
  array[N] int<lower = -B, upper = B> V;        // reported voter positions
  int<lower=0, upper=1> CV;               // indicator of cross-validation
  array[N_obs] int<lower=0, upper=1> holdout;   // holdout for cross-validation
}

transformed data {
  real<lower = 0> sigma_alpha_prior_rate = (2 - 1) / (B / 5.0);
  real<lower = 0> tau_prior_rate = (2 - 1) / (B / 5.0);
  vector<lower = -B, upper = B>[N] Vvec = to_vector(V);
  array[2, N_obs] real<lower = -B, upper = B> p;
  for (n in 1:N_obs) {
    p[1, n] = U[n] * V[ii[n]] + B * U[n] - B;
    p[2, n] = U[n] * V[ii[n]] - B * U[n] + B;
  }
}

parameters {
  matrix[N, 2] alpha_raw;                 // shift parameter, split, raw
  matrix[N, 2] beta_raw;                  // stretch parameter, split, raw
  ordered[2] theta_lr;                    // left and right pole
  array[J] real theta_raw;                      // remaining stimuli
  real<lower = 0> sigma_alpha;            // sd of alpha
  real<lower = 0, upper = 2> sigma_beta;  // sd of log(beta)
  real<lower = 0> tau;                    // sd of errors
  vector<lower = 0, upper = 1>[N] lambda; // mixing proportion, flipping
  real<lower = .5, upper = 1> psi;        // mean of prior on lambda
  real<lower = 2, upper = 100> delta;     // concentration of prior on lambda
  array[N] real<lower = 0, upper = 1> gamma;    // rationalization per respondent
  real<lower = 1> gam_a;                  // hyperparameter for gamma
  real<lower = 1> gam_b;                  // hyperparameter for gamma
  vector<lower = 0, upper = 1>[2 * B + 1] zeta; // direction of rationalization
}

transformed parameters {
  vector[4] log_probs;
  real<lower=0> alpha_lambda = delta * psi; // reparameterization
  real<lower=0> beta_lambda = delta * (1 - psi);
  array[J] real theta;                          // latent stimuli position
  matrix[N, 2] alpha0;                    // shift parameter, split
  matrix[N, 2] beta0;                     // stretch parameter, split
  matrix[N, 2] chi0;                      // latent respondent positions, split
  vector[2] mu0;                          // dif-adjusted mean
  vector[N_obs] log_lik;                  // pointwise log-likelihood for Y
  theta = theta_raw;
  theta[L] = theta_lr[1];                 // safeguard to ensure identification
  theta[R] = theta_lr[2];
  alpha0[, 1] = alpha_raw[, 1] * sigma_alpha; // non-centered specifications
  alpha0[, 2] = alpha_raw[, 2] * sigma_alpha;
  beta0[, 1] = exp(beta_raw[, 1] * sigma_beta);
  beta0[, 2] = -exp(beta_raw[, 2] * sigma_beta);
  chi0[, 1] = ((Vvec - alpha0[, 1]) ./ beta0[, 1]);
  chi0[, 2] = ((Vvec - alpha0[, 2]) ./ beta0[, 2]);

  for (n in 1:N_obs) {
    mu0[1] = alpha0[ii[n], 1] + beta0[ii[n], 1] * theta[jj[n]];
    mu0[2] = alpha0[ii[n], 2] + beta0[ii[n], 2] * theta[jj[n]];
    log_probs[1] = log(lambda[ii[n]]) + log(zeta[V[ii[n]] + B + 1]) +
      normal_lpdf(Y[n] | (1 - gamma[ii[n]]) * mu0[1] + gamma[ii[n]] * p[1, n], tau);
    log_probs[2] = log(lambda[ii[n]]) + log((1 - zeta[V[ii[n]] + B + 1])) +
      normal_lpdf(Y[n] | (1 - gamma[ii[n]]) * mu0[1] + gamma[ii[n]] * p[2, n], tau);
    log_probs[3] = log((1 - lambda[ii[n]])) + log(zeta[V[ii[n]] + B + 1]) +
      normal_lpdf(Y[n] | (1 - gamma[ii[n]]) * mu0[2] + gamma[ii[n]] * p[1, n], tau);
    log_probs[4] = log((1 - lambda[ii[n]])) + log((1 - zeta[V[ii[n]] + B + 1])) +
      normal_lpdf(Y[n] | (1 - gamma[ii[n]]) * mu0[2] + gamma[ii[n]] * p[2, n], tau);
    log_lik[n] = log_sum_exp(log_probs);
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
  gamma ~ beta(gam_a, gam_b);
  gam_a ~ gamma(1.5, .5);
  gam_b ~ gamma(1.5, .5);
  zeta ~ beta(1.2, 1.2);
  tau ~ gamma(2, tau_prior_rate);
  lambda ~ beta(alpha_lambda, beta_lambda);
  psi ~ beta(8.5, 1.5);
  delta - 2 ~ gamma(2, .1);

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
