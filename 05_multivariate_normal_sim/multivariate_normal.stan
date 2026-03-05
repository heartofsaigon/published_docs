data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> D;                       // D=2
  array[N] int<lower=1, upper=K> arm;
  matrix[N, D] Y;

  real<lower=1> eta_y;                  // LKJ shape for outcome correlation
  real<lower=0> mu_scale;               // prior scale for means
  real<lower=0> sd_y_scale;             // prior scale for outcome SDs
}
parameters {
  matrix[K, D] Mu;                      // Mu[k, ] = mean vector for arm k (independent across k a priori)
  vector<lower=1e-8>[D] sd_y;              // marginal SDs for the outcomes
  cholesky_factor_corr[D] L_Omega_y;    // Cholesky factor of outcome correlation
}
transformed parameters {
  matrix[D, D] L_Sigma;
  L_Sigma = diag_pre_multiply(sd_y, L_Omega_y);
}
model {
  // Priors: independent across arms (no hierarchical coupling)
  to_vector(Mu) ~ normal(0, mu_scale);

  // Outcome covariance prior
  sd_y ~ normal(0, sd_y_scale);                 // half-normal via <lower=0>
  L_Omega_y ~ lkj_corr_cholesky(eta_y);

  // Likelihood
  for (i in 1:N) {
    Y[i] ~ multi_normal_cholesky( (Mu[arm[i]])', L_Sigma );
  }
}
generated quantities {
  // posterior predictive replicated data
  matrix[N, D] Y_rep;
  for (i in 1:N) {
    // Draw y_rep[i, ] ~ MVN( Mu[arm[i]], Sigma )
    Y_rep[i] = (multi_normal_cholesky_rng( (Mu[arm[i]])', L_Sigma ))';
  }

}
