
rm(list = ls())
pacman::p_load(tidyverse, MASS, cmdstanr, posterior, HDInterval, magrittr, data.table,
               emoa)

set.seed(44)

# Dimensions
K <- 10         # treatments/arms
D <- 2          # objectives/outcomes per subject
n <- 50        # subjects per arm
N <- K * n      # total subjects



#### helper functions ########

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

#---------------

vikor_score_min <- function(Y, weights = NULL, v = 0.5) {
        # Y: matrix/data.frame, rows = design points, cols = objectives
        # Assumption:
        # - all objectives are cost criteria (smaller is better)
        
        Y <- as.matrix(Y)
        
        if (!is.numeric(Y)) {
                stop("Y must be a numeric matrix or data.frame.")
        }
        
        n <- nrow(Y)
        m <- ncol(Y)
        
        if (is.null(weights)) {
                weights <- rep(1 / m, m)
        }
        
        if (length(weights) != m) {
                stop("Length of weights must equal number of objectives.")
        }
        
        if (any(weights < 0)) {
                stop("weights must be nonnegative.")
        }
        
        if (abs(sum(weights) - 1) > 1e-8) {
                weights <- weights / sum(weights)
        }
        
        if (v < 0 || v > 1) {
                stop("v must be between 0 and 1.")
        }
        
        # For minimization:
        f_star <- apply(Y, 2, min)   # best
        f_minus <- apply(Y, 2, max)  # worst
        
        gap <- matrix(0, nrow = n, ncol = m)
        
        for (j in seq_len(m)) {
                denom <- f_minus[j] - f_star[j]
                
                if (denom == 0) {
                        gap[, j] <- 0
                } else {
                        gap[, j] <- weights[j] * (Y[, j] - f_star[j]) / denom
                }
        }
        
        S <- rowSums(gap)
        R <- apply(gap, 1, max)
        
        S_star <- min(S)
        S_minus <- max(S)
        R_star <- min(R)
        R_minus <- max(R)
        
        if (S_minus == S_star) {
                S_term <- rep(0, n)
        } else {
                S_term <- (S - S_star) / (S_minus - S_star)
        }
        
        if (R_minus == R_star) {
                R_term <- rep(0, n)
        } else {
                R_term <- (R - R_star) / (R_minus - R_star)
        }
        
        Q <- v * S_term + (1 - v) * R_term
        
        out <- data.frame(
                design = seq_len(n),
                S = S,
                R = R,
                Q = Q,
                rank_Q = rank(Q, ties.method = "first")
        )
        
        out <- out[order(out$Q), ]
        rownames(out) <- NULL
        
        return(out)
}
################################################################################


        

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
text(Mu_true[,1], Mu_true[,2]-0.1, col = "red", cex = 0.7)
points(Mu_true[compose(which, Negate(is_dominated), t)(Mu_true),],
       col = "red", cex = .5, pch = 16)
# fitting model ----------------------------


mod <- cmdstan_model(stan_file = "/Users/nam-anhtran/heartofsaigon/published_docs/05_multivariate_normal_sim/multivariate_normal.stan")

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
        rbind(i[1:K], i[(K+1):(K*2)])%>%
                Negate(is_dominated)()%>%
                which()
}, simplify = FALSE)

# probability of being nondominated using posteriors
sapply(1:10, \(i) sapply(dom, \(k) i%in%k)|> mean())|>
        `names<-`(1:10)%T>%
        {barplot(.)}


########

idx_non<-
rbind(est[1:K,]$mean, est[(K+1):(2*K),]$mean)|>
        Negate(is_dominated)()|>
        which()


r<- 
sapply(1:nrow(sam_y1), \(i){
        d = rbind(sam_y1[i,idx_non], sam_y2[i,idx_non])%>% t()
        `rownames<-`(d, c(1:4,7))
        vikor_score_min(d)[,"design"]
})|> t()

r[,1]%>% table()

