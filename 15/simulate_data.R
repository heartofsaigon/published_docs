

library(cmdstanr)
library(posterior)
library(magrittr)
library(bayesplot)

set.seed(20260331)

## ---------------------------
## 1. Settings
## ---------------------------
K <- 10
m <- 50
D <- 2
N <- K * m

## ---------------------------
## 2. True treatment means
##    Each row is one treatment,
##    each column is one outcome
## ---------------------------
mu_true <- -matrix(
  c(
    # --- Pareto (nondominated) points ---
    0.5, 2.5,
    0.9, 2.2,
    1.3, 1.8,
    1.8, 1.3,
    2.3, 0.9,
    
    # --- Dominated points ---
    0.8, 1.8,   # dominated by (0.9, 2.2)
    1.2, 1.2,   # dominated by (1.3, 1.8)
    1.6, 1.0,   # dominated by (1.8, 1.3)
    1.0, 1.5,   # dominated by (1.3, 1.8)
    2.0, 0.8    # dominated by (2.3, 0.9)
  ),
  nrow = K,
  ncol = D,
  byrow = TRUE
)

## ---------------------------
## 3. True covariance matrix
## ---------------------------
sigma_true <- c(0.8, 1.1)
rho_true <- 0.45

Sigma_true <- matrix(
 c(
  sigma_true[1]^2,
  rho_true * sigma_true[1] * sigma_true[2],
  rho_true * sigma_true[1] * sigma_true[2],
  sigma_true[2]^2
 ),
 nrow = D,
 ncol = D,
 byrow = TRUE
)

L_true <- t(chol(Sigma_true))

## ---------------------------
## 4. Treatment assignment
## ---------------------------
trt <- rep(1:K, each = m)

## ---------------------------
## 5. Simulate outcomes
##    y_i = mu_{trt_i} + e_i
##    e_i ~ N_2(0, Sigma_true)
## ---------------------------
Y_mat <- matrix(NA_real_, nrow = N, ncol = D)

for (i in 1:N) {
 eps_i <- drop(L_true %*% rnorm(D))
 Y_mat[i, ] <- mu_true[trt[i], ] + eps_i
}

colnames(Y_mat) <- c("y1", "y2")

## ---------------------------
## 7. Create Stan data list
## ---------------------------
stan_data <- list(
 N = N,
 K = K,
 D = D,
 trt = trt,
 y = Y_mat
)

## ---------------------------
## 8. Compile the Stan model
## ---------------------------
mod <- cmdstan_model("15/multinormal.stan")

## ---------------------------
## 9. Fit the model
## ---------------------------
fit <- mod$sample(
 data = stan_data,
 seed = 20260331,
 chains = 4,
 parallel_chains = 4,
 iter_warmup = 1500,
 iter_sampling = 1000,
 refresh = 200
)

## ---------------------------
## 10. Basic summaries
## ---------------------------

est_mu<-
fit$summary(c("mu"))|>
  dplyr::pull(mean)|>
  matrix(ncol = 2, nrow = K, byrow = FALSE)

purrr::compose(which, Negate(emoa::is_dominated), t)(est_mu)
purrr::compose(which, Negate(emoa::is_dominated), t)(mu_true)
## ---------------------------
## 11. Extract posterior means
## ---------------------------
draws_df <- fit$draws("mu", format = "matrix")%>%
  apply(1, function(x){
    m = matrix(x, ncol = D, nrow = K) 
    purrr::compose(which, Negate(emoa::is_dominated), t)(m)
    }, simplify = F)

sapply(1:10, \(k) mean(sapply(draws_df, \(K) k%in%K)))


fit$summary("mu")$mean%>% matrix(ncol = D, nrow = K)
mu_true

#### diagnostic check 

draws = fit$draws("mu", format = "matrix")
mcmc_trace(draws)
mcmc_dens_overlay(draws)
mcmc_rank_overlay(draws)


colnames(draws)





