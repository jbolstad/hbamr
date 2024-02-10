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
  vector[J] mean_spos;                    // average stimuli placements
  int<lower = 0, upper = 1> CV;           // indicator of cross-validation
  vector<lower = 0, upper = 1>[N_obs] holdout; // holdout for cross-validation
}

transformed data {
  real<lower = 0> tau_prior_rate = (2 - 1) / (B / 5.0);
  vector<lower = 0, upper = 1>[N_obs] not_holdout = 1 - holdout;
}

parameters {
  vector[N] alpha;                        // shift parameter
  vector[N] beta;                         // stretch parameter
  real<lower = mean_spos[L] - B / 50.0, upper = mean_spos[L]> thetal; // left pole stimuli
  real<lower = mean_spos[R], upper = mean_spos[R] + B / 50.0> thetar; // right pole stimuli
  array[J] real theta_raw;                // remaining stimuli
  real<lower = 3, upper = 30> nu;         // concentration of etas
  real<lower = 0> tau;                    // scale of errors
  vector<lower = 0>[N] eta;               // mean ind. error variance x J^2
  simplex[J] rho;                         // stimuli-shares of variance
}

transformed parameters {
  array[J] real theta;                    // latent stimuli position
  vector[N_obs] log_lik;                  // pointwise log-likelihood
  real<lower = 0> eta_scale = tau * J;
  theta = theta_raw;
  theta[L] = thetal;
  theta[R] = thetar;
  for (n in 1:N_obs) {
    log_lik[n] = normal_lpdf(Y[n] | alpha[ii[n]] + beta[ii[n]] * theta[jj[n]],
      sqrt(eta[ii[n]]) * rho[jj[n]]);
  }
}

model {
  alpha ~ uniform(-100, 100);
  beta ~ uniform(-100, 100);
  theta_raw ~ normal(0, B);
  thetal ~ normal(0, B);
  thetar ~ normal(0, B);
  eta ~ scaled_inv_chi_square(nu, eta_scale);
  nu ~ gamma(25, 2.5);
  tau ~ gamma(2, 5 / (B * 1.0));
  rho ~ dirichlet(rep_vector(5, J));

  if (CV == 0)
    target += sum(log_lik);
  else
    target += sum(log_lik .* not_holdout);
}

generated quantities {
  real<lower = 0> min_rho = min(rho);
  vector[N] chi = ((V - to_vector(normal_rng(0, sqrt(eta) * min_rho)) - alpha) ./ beta);
}
