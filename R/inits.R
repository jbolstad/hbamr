inits_HBAM <- function(chain_id = 1, dat) {
  list (
    theta_raw = dat$mean_spos + rnorm(dat$J, 0, (dat$B / 5) * 0.25),
    theta_lr = dat$mean_spos[c(dat$L, dat$R)] + c(runif(1, -dat$B / 10, 0), runif(1, 0, dat$B / 10)),
    sigma_chi = rinvchisq(1, 500, dat$B / 2),
    chi0 = matrix(rnorm(dat$N * 2, 0, dat$B / 4), ncol = 2),
    sigma_alpha = rinvchisq(1, 100, dat$B / 5),
    alpha_raw = matrix(rnorm(2 * dat$N, 0, .5), ncol = 2),
    sigma_beta = runif(1, .1, .2),
    beta_raw = matrix(rnorm(dat$N * 2, 0, .2), ncol = 2),
    lambda = runif(dat$N, .7, .95),
    nu = 3 + rinvchisq(1, 100, 7),
    tau = rinvchisq(1, 500, (dat$B / 3)),
    eta = rinvchisq(dat$N, 100, dat$J^2 * (dat$B / 3)^2),
    rho = c(rdirichlet(1, rep(50, dat$J))),
    logit_lambda = rnorm(dat$N, 0, .2),
    psi = exp(rnorm(1, 1.4, .1)),
    # For HBAM_ORIG:
    delta = runif(1, 2.05, 2.1),
    psi_old = runif(1, .875, .925),
    # For HBAM_MULTI:
    mu_alpha_raw = rdirichlet(1, rep(500, dat$G)),
    mu_beta_raw = rdirichlet(1, rep(500, dat$G)),
    # For HBAM_R:
    gammma = runif(dat$N, 0, .3),
    gam_a = runif(1, 2, 3),
    gam_b = runif(1, 1.1, 1.3),
    zeta = runif(dat$B * 2 + 1, .1, .9)
  )
}

inits_HBAM_NF <- function(chain_id = 1, dat) {
  list (
    theta_raw = dat$mean_spos + rnorm(dat$J, 0, (dat$B / 5) * 0.25),
    theta_lr = dat$mean_spos[c(dat$L, dat$R)] + c(runif(1, -dat$B / 10, 0), runif(1, 0, dat$B / 10)),
    sigma_chi = rinvchisq(1, 500, dat$B / 2),
    chi = rnorm(dat$N, 0, dat$B / 4),
    sigma_alpha = rinvchisq(1, 100, dat$B / 5),
    alpha_raw = rnorm(dat$N, 0, .5),
    sigma_beta = runif(1, .1, .2),
    beta_raw = rnorm(dat$N, 0, .2),
    nu = 3 + rinvchisq(1, 100, 7),
    tau = rinvchisq(1, 500, (dat$B / 3)),
    eta = rinvchisq(dat$N, 100, dat$J^2 * (dat$B / 3)^2),
    rho = c(rdirichlet(1, rep(50, dat$J))),
    # For HBAM_MULTI:
    mu_alpha_raw = rdirichlet(1, rep(500, dat$G)),
    mu_beta_raw = rdirichlet(1, rep(500, dat$G))
  )
}

inits_BAM <- function(chain_id = 1, dat) {
  list (
    theta_raw = dat$mean_spos + rnorm(dat$J, 0, (dat$B / 5) * 0.25),
    beta = rnorm(dat$N, 1, .5),
    thetal = -1.05, thetar = 1.05,
    nu = 3 + rinvchisq(1, 100, 7),
    tau = rinvchisq(1, 500, (dat$B / 3)),
    eta = rinvchisq(dat$N, 100, dat$J^2 * (dat$B / 3)^2),
    rho = c(rdirichlet(1, rep(50, dat$J)))
  )
}

# Collecting all inits-functions:
inits <- list(
  HBAM = inits_HBAM,
  HBAM_ORIG = inits_HBAM,
  HBAM_MAX = inits_HBAM,
  # HBAM_2 = inits_HBAM,
  # HBAM_HM = inits_HBAM,
  HBAM_MINI = inits_HBAM,
  HBAM_MULTI = inits_HBAM,
  # HBAM_R = inits_HBAM,
  HBAM_R_MINI = inits_HBAM,
  FBAM_MINI = inits_HBAM,
  FBAM_MULTI = inits_HBAM,
  HBAM_NF = inits_HBAM_NF,
  HBAM_MULTI_NF = inits_HBAM_NF,
  FBAM_MULTI_NF = inits_HBAM_NF,
  BAM = inits_BAM)

rinvchisq <- function(n, df, scale = 1/df) {
  return((df * scale) / rchisq(n, df = df))
}

rdirichlet <- function (n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x / as.vector(sm))
}
