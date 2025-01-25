data {
  int<lower = 1> N;                       // n of individuals
  int<lower = 1> J;                       // n of items
  int<lower = 0> G;                       // n of groups
  int<lower = 1> N_obs;                   // n of observations
  array[N_obs] int<lower = 1> ii;         // index i in matrix
  array[N_obs] int<lower = 1> jj;         // index j in matrix
  array[N] int<lower = 1> gg;             // group-index of i
  int<lower = 1> B;                       // length of scale -1 / 2
  int<lower = 1, upper = J> L;            // left pole
  int<lower = 1, upper = J> R;            // right pole
  array[N_obs] int<lower = -B, upper = B> Y; // reported stimuli positions
  array[N_obs] real<lower = 0, upper = 1> U; // reported voter preferences
  array[N] int<lower = -B, upper = B> V;  // reported voter positions
  vector[J] mean_spos;                    // average stimuli placements
  real<lower = 0> sigma_mu_alpha;         // sd of prior on mu_alpha, MULTI
  real<lower = 0> sigma_mu_beta;          // sd of prior on mu_beta, MULTI
  real<lower = 0> sigma_alpha_fixed;      // sd of prior on alpha, FBAM
  real<lower = 0> sigma_beta_fixed;       // sd of prior on log(beta), FBAM
  int<lower = 0, upper = 1> flip;         // allow flipping?
  int<lower = 0, upper = 1> het;          // model heteroskedasticity?
  int<lower = 0, upper = 1> rat;          // model rationalization?
  int<lower = 0, upper = 1> group;        // give groups different hyperpriors?
  int<lower = 0, upper = 1> bam;          // unpooled model with uniform priors?
  int<lower = 0, upper = 1> fixed;        // use fixed priors, not hyperparams.
  int<lower = 0, upper = 1> CV;           // indicator of cross-validation
  int<lower = 0, upper = 1> MCMC;         // indicator of fitting method
  vector<lower = 0, upper = 1>[N_obs] holdout; // holdout for cross-validation
}

transformed data {
  real nu_fixed = 6;                      // concentration of etas
  real tau_fixed = B / 4.0;               // scale of prior on errors
  real eta_scale_fixed = tau_fixed * J;
  real psi_fixed = 6;                     // implies 13% prior prob. of flipping
  real<lower = 0> sigma_alpha_prior_rate = (5 - 1) / (B / 8.0);
  real<lower = 0> tau_prior_rate = (2 - 1) / (B / 5.0);
  vector<lower = 0, upper = 1>[N_obs] not_holdout = 1 - holdout;
  real mean_mu_simplexes = 1.0 / G;       // for later scaling of simplexes
  real sd_mu_simplexes = sqrt(mean_mu_simplexes * (1 - mean_mu_simplexes) / (50 * G + 1));
  array[2, N_obs] real<lower = -B, upper = B> p;
  for (n in 1:N_obs) {
    p[1, n] = U[n] * V[ii[n]] + B * U[n] - B;
    p[2, n] = U[n] * V[ii[n]] - B * U[n] + B;
  }
}

parameters {
  matrix[N, 1 + flip] alpha_raw;          // shift parameter, split, raw
  matrix[N, 1 + flip] beta_raw;           // stretch parameter, split, raw
  ordered[bam == 0 ? 2 : 0] theta_lr;     // left and right pole
  array[J] real theta_raw;                // remaining stimuli
  simplex[group == 1 ? G : 1] mu_alpha_raw; // group-level mean of alpha, raw
  simplex[group == 1 ? G : 1] mu_beta_raw; // group-level mean of log(beta), raw
  array[bam == 0 && fixed == 0 ? 1 : 0] real<lower = 0> sigma_alpha_par; // sd of alpha
  array[bam == 0 && fixed == 0 ? 1 : 0] real<lower = 0, upper = 2> sigma_beta_par; // sd of log(beta)
  array[het == 1 && fixed == 0 ? 1 : 0] real<lower = 3, upper = 30> nu_par; // concentration of etas
  array[1 - fixed] real<lower = 0> tau_par; // scale of errors
  vector<lower = 0>[het == 1 ? N : 0] eta; // mean ind. error variance x J^2
  simplex[het == 1 ? J : 1] rho;          // stimuli-shares of variance
  vector[flip == 1 ? N : 0] lambda_raw;   // raw mixing proportion, flipping
  array[flip == 1 && fixed == 0 ? 1 : 0] real<lower = 0> psi_par; // mean of prior on logit of lambda
  array[rat == 1 ? N : 0] real<lower = 0, upper = 1> gamma; // rationalization per respondent
  array[rat] real<lower = 1> gam_a;       // hyperparameter for gamma
  array[rat] real<lower = 1> gam_b;       // hyperparameter for gamma
  vector<lower = 0, upper = 1>[rat == 1 ? 2 * B + 1 : 0] zeta; // direction of rationalization
}

transformed parameters {
  vector[rat == 1 ? 4 : 0] log_probs;
  array[J] real theta;                    // latent stimuli position
  matrix[N, 1 + flip] alpha0;             // shift parameter, split
  matrix[N, 1 + flip] beta0;              // stretch parameter, split
  vector[N_obs] log_lik;                  // pointwise log-likelihood for Y
  vector<lower = 0, upper = 1>[flip == 1 ? N : 0] lambda;
  vector[group == 1 ? G : 0] mu_alpha;
  vector[group == 1 ? G : 0] mu_beta;
  array[het] real<lower = 0> eta_scale;
  array[1 - bam] real<lower = 0> sigma_alpha; // sd of alpha
  array[1 - bam] real<lower = 0, upper = 2> sigma_beta; // sd of log(beta)
  array[het] real<lower = 3, upper = 30> nu; // concentration of etas
  real<lower = 0> tau;                    // scale of errors
  array[flip] real<lower = 0> psi;        // mean of prior on logit of lambda

  if (fixed == 0) {
    tau = tau_par[1];
    if (bam == 0) {
      sigma_alpha = sigma_alpha_par;
      sigma_beta = sigma_beta_par;
    }
    if (het == 1) {
      nu = nu_par;
      eta_scale[1] = tau * J;
    }
    if (flip == 1) {
      psi = psi_par;
    }
  } else {
    tau = tau_fixed;
    if (bam == 0) {
      sigma_alpha[1] = sigma_alpha_fixed;
      sigma_beta[1] = sigma_beta_fixed;
    }
    if (het == 1) {
      nu[1] = nu_fixed;
      eta_scale[1] = tau * J;
    }
    if (flip == 1) {
      psi[1] = psi_fixed;
    }
  }
  if (group == 1) {
    mu_alpha = ((mu_alpha_raw - mean_mu_simplexes) / sd_mu_simplexes) * sigma_mu_alpha;
    mu_beta = ((mu_beta_raw - mean_mu_simplexes) / sd_mu_simplexes) * sigma_mu_beta;
  }
  theta = theta_raw;

  if (bam == 1) {
    theta[L] = mean_spos[L];
    theta[R] = mean_spos[R];
    alpha0[, 1] = alpha_raw[, 1];
    beta0[, 1] = beta_raw[, 1];
    for (n in 1:N_obs) {
      log_lik[n] = normal_lpdf(Y[n] | alpha0[ii[n], 1] + beta0[ii[n], 1] * theta[jj[n]],
        sqrt(eta[ii[n]]) * rho[jj[n]]);
    }
  } else {
    theta[L] = theta_lr[1];               // safeguard to ensure identification
    theta[R] = theta_lr[2];
    if (group == 0) {
      alpha0[, 1] = alpha_raw[, 1] * sigma_alpha[1];
      beta0[, 1] = exp(beta_raw[, 1] * sigma_beta[1]);
    } else {
      alpha0[, 1] = alpha_raw[, 1] * sigma_alpha[1] + mu_alpha[gg];
      beta0[, 1] = exp(beta_raw[, 1] * sigma_beta[1] + mu_beta[gg]);
    }
    if (flip == 1) {
      lambda = inv_logit(psi[1] + lambda_raw * 3); // prob. of non-flipping
      beta0[, 2] = -exp(beta_raw[, 2] * sigma_beta[1]);
      if (group == 0) {
        alpha0[, 2] = alpha_raw[, 2] * sigma_alpha[1];
      } else {
        alpha0[, 2] = alpha_raw[, 2] * sigma_alpha[1] + mu_alpha[gg];
      }
    }
    if (rat == 1) {                       // HBAM_R_MINI
      for (n in 1:N_obs) {
        vector[2] mu0;                    // dif-adjusted mean
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
    } else if (het == 0) {                // HBAM_MINI
      for (n in 1:N_obs) {
        log_lik[n] = log_mix( lambda[ii[n]],
          normal_lpdf(Y[n] | alpha0[ii[n], 1] + beta0[ii[n], 1] * theta[jj[n]], tau),
          normal_lpdf(Y[n] | alpha0[ii[n], 2] + beta0[ii[n], 2] * theta[jj[n]], tau) );
      }
    } else {
      if (flip == 1) {                    // HBAM & HBAM_MULTI
        for (n in 1:N_obs) {
          log_lik[n] = log_mix( lambda[ii[n]],
            normal_lpdf(Y[n] | alpha0[ii[n], 1] + beta0[ii[n], 1] * theta[jj[n]],
              sqrt(eta[ii[n]]) * rho[jj[n]]),
            normal_lpdf(Y[n] | alpha0[ii[n], 2] + beta0[ii[n], 2] * theta[jj[n]],
              sqrt(eta[ii[n]]) * rho[jj[n]]) );
        }
      } else {                            // HBAM_NF & HBAM_MULTI_NF
        for (n in 1:N_obs) {
          log_lik[n] = normal_lpdf(Y[n] | alpha0[ii[n], 1] + beta0[ii[n], 1] * theta[jj[n]],
            sqrt(eta[ii[n]]) * rho[jj[n]]);
        }
      }
    }
  }
}

model {
  theta_raw ~ normal(0, B / 2.0);
  if (bam == 1) {
    alpha_raw[, 1] ~ uniform(-100, 100);
    beta_raw[, 1] ~ uniform(-100, 100);
  } else {
    theta_lr ~ normal(0, B / 2.0);
    alpha_raw[, 1] ~ normal(0, 1);
    sigma_alpha ~ gamma(5, sigma_alpha_prior_rate);
    beta_raw[, 1] ~ normal(0, 1);
    sigma_beta ~ gamma(9, 40);
  }
  if (group == 1) {
    mu_alpha_raw ~ dirichlet(rep_vector(50, G));
    mu_beta_raw ~ dirichlet(rep_vector(50, G));
  } else {
    mu_alpha_raw ~ dirichlet(rep_vector(50, 1));
    mu_beta_raw ~ dirichlet(rep_vector(50, 1));
  }
  if (flip == 1) {
    alpha_raw[, 2] ~ normal(0, 1);
    beta_raw[, 2] ~ normal(0, 1);
    lambda_raw ~ normal(0, 1);
    psi ~ lognormal(1.4, .5);
  }
  if (rat == 1) {
    gamma ~ beta(gam_a[1], gam_b[1]);
    gam_a ~ gamma(1.5, .5);
    gam_b ~ gamma(1.5, .5);
    zeta ~ beta(1.2, 1.2);
  }
  if (het == 1) {
    eta ~ scaled_inv_chi_square(nu[1], eta_scale[1]);
    nu ~ gamma(25, 2.5);
    rho ~ dirichlet(rep_vector(50, J));
  } else {
    rho ~ dirichlet(rep_vector(50, 1));
  }
  tau ~ gamma(2, tau_prior_rate);

  if (CV == 0)
    target += sum(log_lik);
  else
    target += sum(log_lik .* not_holdout);
}

generated quantities {
  real<lower = 0> min_rho = min(rho);
  vector[flip == 1 ? N : 0] kappa;
  vector[N] alpha;
  vector[N] beta;
  vector[N] chi;
  vector[rat == 0 && MCMC == 1 ? N_obs : 0] Y_pred;
  if (flip == 1) {
    if (MCMC == 1)
      kappa = to_vector(bernoulli_rng(lambda));
    else
      kappa = to_vector(round(lambda));   // Rounding to MAP instead of sampling
    alpha = (kappa .* alpha0[, 1]) + ((1 - kappa) .* alpha0[, 2]);
    beta = (kappa .* beta0[, 1]) + ((1 - kappa) .* beta0[, 2]);
  } else {
    alpha = alpha0[, 1];
    beta = beta0[, 1];
  }
  if (het == 1) {
    if (MCMC == 1) {
      chi = (to_vector(V) - to_vector(normal_rng(0, sqrt(eta) * min_rho)) - alpha) ./ beta;
      for (n in 1:N_obs)
        Y_pred[n] = normal_rng(alpha[ii[n]] + beta[ii[n]] * theta[jj[n]], sqrt(eta[ii[n]]) * rho[jj[n]]);
    }
    else {
      chi = (to_vector(V) - alpha) ./ beta;
    }
  } else {
    if (MCMC == 1) {
      chi = (to_vector(V) - to_vector(normal_rng(0, rep_vector(tau, N))) - alpha) ./ beta;
      if (rat == 0) {
        for (n in 1:N_obs)
          Y_pred[n] = normal_rng(alpha[ii[n]] + beta[ii[n]] * theta[jj[n]], tau);
      }
    } else {
      chi = (to_vector(V) - alpha) ./ beta;
    }
  }
}
