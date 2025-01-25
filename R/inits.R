inits_omni <- function(chain_id = 1, dat) {
  alpha_raw <- matrix(rnorm(2 * dat$N, 0, .3), ncol = 2)
  beta_raw <- matrix(rnorm(dat$N * 2, 0, .2), ncol = 2)
  if (dat$fixed == 0) {
    tau <- rinvchisq(1, 500, (dat$B / 3))
    dim(tau) <- 1
  } else {
    tau <- NULL
  }
  if (dat$bam == 0 & dat$fixed == 0) {
    theta_lr <- dat$mean_spos[c(dat$L, dat$R)] + c(runif(1, -dat$B / 10, 0), runif(1, 0, dat$B / 10))
    sigma_alpha <- rinvchisq(1, 200, .75)
    sigma_beta <- runif(1, .3, .35)
    dim(sigma_alpha) <- 1
    dim(sigma_beta) <- 1
  } else {
    theta_lr <- NULL
    sigma_alpha <- NULL
    sigma_beta <- NULL
  }
  if (dat$fixed == 1) {
    theta_lr <- dat$mean_spos[c(dat$L, dat$R)] + c(runif(1, -dat$B / 10, 0), runif(1, 0, dat$B / 10))
  }
  if (dat$flip == 0) {
    alpha_raw <- alpha_raw[, 1]
    dim(alpha_raw) <- c(dat$N, 1)
    beta_raw <- beta_raw[, 1]
    dim(beta_raw) <- c(dat$N, 1)
    lambda_raw <- NULL
    psi <- NULL
  } else {
    lambda_raw <- rnorm(dat$N, 0, .05)
    if (dat$fixed == 0) {
      psi <- exp(rnorm(1, 1.3, .1))
      dim(psi) <- 1
    } else {
      psi <- NULL
    }
  }
  if (dat$het == 1) {
    nu <- 3 + rinvchisq(1, 100, 7)
    eta <- rinvchisq(dat$N, 100, dat$J^2 * (dat$B / 3)^2)
    rho <- rdirichlet(1, rep(75, dat$J))
    dim(nu) <- 1
  } else {
    nu <- NULL
    eta <- NULL
    rho <- 1
    dim(rho) <- 1
  }
  if (dat$fixed == 1) {
    nu <- NULL
  }
  if (dat$group == 1) {
    mu_alpha_raw <- rdirichlet(1, rep(500, dat$G))
    mu_beta_raw <- rdirichlet(1, rep(500, dat$G))
  } else {
    mu_alpha_raw <- 1
    mu_beta_raw <- 1
    dim(mu_alpha_raw) <- 1
    dim(mu_beta_raw) <- 1
  }
  if (dat$rat == 1) {
    gammma <- runif(dat$N, 0, .3)
    gam_a <- runif(1, 2, 3)
    gam_b <- runif(1, 1.1, 1.3)
    dim(gam_a) <- 1
    dim(gam_b) <- 1
    zeta <- runif(dat$B * 2 + 1, .1, .9)
  } else {
    gammma <- NULL
    gam_a <- NULL
    gam_b <- NULL
    zeta <- NULL
  }
  list (
    theta_raw = dat$mean_spos + rnorm(dat$J, 0, (dat$B / 5) * 0.25),
    theta_lr = theta_lr,
    sigma_alpha_par = sigma_alpha,
    alpha_raw = alpha_raw,
    sigma_beta_par = sigma_beta,
    beta_raw = beta_raw,
    nu_par = nu,
    tau_par = tau,
    eta = eta,
    rho = rho,
    lambda_raw = lambda_raw,
    psi_par = psi,
    mu_alpha_raw = mu_alpha_raw,
    mu_beta_raw = mu_beta_raw,
    gammma = gammma,
    gam_a = gam_a,
    gam_b = gam_b,
    zeta = zeta
  )
}

rinvchisq <- function(n, df, scale = 1/df) {
  return((df * scale) / rchisq(n, df = df))
}

rdirichlet <- function (n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  out <- x / as.vector(sm)
  if (n == 1) { out <- as.numeric(out) }
  return(out)
}
