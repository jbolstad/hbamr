# Function to draw from the scaled inverse chi-squared distribution:
rinvchisq <- function(n, df, scale = 1 / df) {
  (df * scale) / rchisq(n, df = df)
}

# Function to draw from the Dirichlet distribution:
rdirichlet <- function (n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x / as.vector(sm))
}

# Function to generate data:
create_data <- function(N = 500, J = 4, B = 3, nu = 8.4, tau = .19 * B, rho_a = 5, psi = .89,
                        sigma_a = .19 * B, sigma_b = .27, sigma_c = .52 * B,
                        seed = 1) {
  set.seed(seed)

  # Latent parameters:
  true_alpha <- rnorm(N, 0, sigma_a)
  true_beta0 <- exp(rnorm(N, 0, sigma_b))
  true_kappa <- rbinom(N, 1, psi)
  true_beta <- true_beta0 * true_kappa - true_beta0 * (1 - true_kappa)
  true_chi <- rnorm(N, 0, sigma_c)
  true_theta <- c(runif(1, -B + 1, 0), runif(1, 0, B - 1), runif(J - 2, -B + 1, B - 1))

  # Stimulus positions:
  true_rho <- rdirichlet(1, rep(rho_a, J))
  true_eta <- rinvchisq(N, nu, scale = (tau^2) * (J^2))
  sd_ij <- matrix(rep(sqrt(true_eta), times = J) * rep(true_rho, each = N), nrow = N, ncol = J)
  errors <- matrix(rnorm(N * J, 0, sd_ij), nrow = N, ncol = J)
  ppos <- round(matrix(rep(true_alpha, each = J) + (rep(true_beta, each = J) * true_theta),
                       N, J, byrow = TRUE) + errors)
  ppos[ppos < -B] <- -B; ppos[ppos > B] <- B

  # Self-placements:
  self <- round(true_alpha + true_beta * true_chi + rnorm(N, 0, sqrt(true_eta) * min(true_rho)))
  self[self < -B] <- -B; self[self > B] <- B

  dat <- prep_data(self, ppos, B = B)
  dat$true_theta <- true_theta
  dat$true_chi <- true_chi[dat$keep]
  dat$true_alpha <- true_alpha[dat$keep]
  dat$true_beta <- true_beta[dat$keep]
  dat$true_kappa <- true_kappa[dat$keep]

  return(dat)
}
