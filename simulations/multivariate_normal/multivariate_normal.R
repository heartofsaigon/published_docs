
rm(list = ls())
pacman::p_load(tidyverse, MASS, cmdstanr, posterior, HDInterval, magrittr, data.table,
               emoa)

set.seed(4444)

# Dimensions
K <- 5         # treatments/arms
D <- 2          # objectives/outcomes per subject
n <- 50        # subjects per arm
N <- K * n      # total subjects

# Arm-specific mean vectors mu_k (K x D) - choose any values
dtlz2_m2 <- function(X) {
        # DTLZ2 with M = 2 objectives (minimization)
        # X: N x d matrix; here you have N = 5 rows, d = 2 columns, entries in [0,1]
        if (!is.matrix(X)) X <- as.matrix(X)
        if (ncol(X) < 2L) stop("Need at least 2 columns for M = 2.")
        
        n <- nrow(X)
        d <- ncol(X)
        
        # g(x2,...,xd) = sum_{i=2}^d (x_i - 0.5)^2
        g <- if (d == 2L) (X[, 2] - 0.5)^2 else rowSums((X[, 2:d, drop = FALSE] - 0.5)^2)
        one_plus_g <- 1 + g
        
        half_pi <- pi / 2
        x1 <- X[, 1]
        
        f1 <- one_plus_g * cos(x1 * half_pi)
        f2 <- one_plus_g * sin(x1 * half_pi)
        
        cbind(f1 = f1, f2 = f2)
}
Mu_true <- -dtlz2_m2(runif(D*K)|> matrix(ncol =D))

t(Mu_true)%>% Negate(is_dominated)()%>% which()

        

# Within-subject covariance across the 2 outcomes (shared across all arms)
rho <- 0.5
sd_y <- c(0.1, 0.3)

Omega <- matrix(rho, D, D); diag(Omega) <- 1
Sigma <- diag(sd_y) %*% Omega %*% diag(sd_y)   # 2x2 covariance

# Simulate: subjects independent, arms independent given Mu_true
arm_id <- rep(1:K, each = n)
Y <- matrix(NA_real_, nrow = N, ncol = D)

for (i in 1:N) {
 k <- arm_id[i]
 Y[i, ] <- MASS::mvrnorm(n = 1, mu = Mu_true[k, ], Sigma = Sigma)
}

colnames(Y) <- c("y1", "y2")
dat <- data.frame(id = 1:N, arm = arm_id, Y)

# Optional: show empirical within-subject correlation pooled over all subjects
print(cor(Y))

plot(dat[,3:4], col = "gray")
points(Mu_true, col = "black")
points(Mu_true[compose(which, Negate(is_dominated), t)(Mu_true),],
       col = "red", cex = .5, pch = 16)
# fitting model ----------------------------


mod <- cmdstan_model(stan_file = "simulations/multivariate_normal/multivariate_normal.stan")

stan_data <- list(
 N = nrow(dat),
 K = K,
 D = D,
 arm = dat$arm,
 Y = as.matrix(dat[, c("y1", "y2")]),
 eta_y = 2.0,
 mu_scale = 5.0,
 sd_y_scale = 1.0
)

fit <- mod$sample(
 data = stan_data,
 chains = 4,
 parallel_chains = 4,
 iter_warmup = 2000,
 iter_sampling = 1000,
 refresh = 200
)

# Posterior summaries for Mu[1..K, 1..D]
est = fit$summary("Mu", compose(~`names<-`(.,"mean"), mean), hdi, rhat);est
est_map = mod$optimize(data = stan_data, jacobian      = TRUE)
sam = fit$draws("Y_rep", format = "matrix")
sam_y1 = sam[,1:(K*n)]; sam_y2 = sam[,(K*n+1):(K*n*2)]
#---


# nondominated using mean
rbind(est$mean[1:K], est$mean[(K+1):(K*2)])%>%
        Negate(is_dominated)()%>%
        which()
est_map$summary("Mu")$estimate%>%
        {rbind(.[1:K], .[(K+1):(K*2)])}%>%
        Negate(is_dominated)()%>%
        which()
# true nondominated
Mu_true%>% t()%>% Negate(is_dominated)()%>% which()


sam_par = fit$draws("Mu", format = "matrix")
dom<- 
apply(sam_par, 1, \(i){
        rbind(i[1:K], i[(K+1):(k*2)])%>%
                Negate(is_dominated)()%>%
                which()
}, simplify = FALSE)

# probability of being nondominated using posteriors
sapply(1:5, \(i) sapply(dom, \(k) i%in%k)|> mean())





