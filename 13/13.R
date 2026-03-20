
pacman::p_load(cmdstanr, HDInterval, tidyverse, posterior)
set.seed(123)

# -----------------------------
# Simulate a typical SEM dataset
# -----------------------------
n <- 200

# Observed covariate
x <- rnorm(n, mean = 0, sd = 1)

# True structural parameters
alpha_true <- 0.5
beta_true  <- 0.8
tau_true   <- 0.7

# Latent depression
eta <- rnorm(n, mean = alpha_true + beta_true * x, sd = tau_true)

# True measurement parameters
nu_true <- c(1.0, 1.5, 0.5)
lambda_true <- c(1.0, 1.2, 0.8)   # lambda_1 fixed to 1 for identification
sigma_true <- c(0.5, 0.6, 0.4)

# Observed indicators
y1 <- rnorm(n, mean = nu_true[1] + lambda_true[1] * eta, sd = sigma_true[1])
y2 <- rnorm(n, mean = nu_true[2] + lambda_true[2] * eta, sd = sigma_true[2])
y3 <- rnorm(n, mean = nu_true[3] + lambda_true[3] * eta, sd = sigma_true[3])

dat <- data.frame(
 id = 1:n,
 x = x,
 y1 = y1,
 y2 = y2,
 y3 = y3
)


glimpse(dat)
# -----------------------------
# Stan model for Bayesian SEM
# -----------------------------

mod = cmdstanr::cmdstan_model(stan_file = "stancode.stan")

stan_data <- list(
 N = nrow(dat),
 x = dat$x,
 y1 = dat$y1,
 y2 = dat$y2,
 y3 = dat$y3
)

fit <- mod$sample(
 data = stan_data,
 seed = 123,
 chains = 4,
 parallel_chains = 4,
 iter_warmup = 2000,
 iter_sampling = 1000,
 refresh = 200)

fit$summary(c('beta', 'tau', 'nu1', 'nu2', 'nu3', 'lambda2', 'lambda3',
              'sigma1', 'sigma2', 'sigma3'), mean, hdi, rhat)
