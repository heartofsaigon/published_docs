data {
  int<lower=1> N;                         // total number of subjects
  int<lower=1> K;                         // number of treatments
  int<lower=1> D;                         // number of outcomes
  array[N] int<lower=1, upper=K> trt;     // treatment index
  matrix[N,D] y;                   // multivariate outcomes
}

parameters {
  vector[D] mu0;                          // population mean across treatments

  array[K] vector[D] z_mu;                // non-centered latent treatment effects

  vector<lower=1e-8>[D] sigma_within;     // within-subject SDs
  cholesky_factor_corr[D] Lcorr_within;   // within-subject correlation

  vector<lower=1e-8>[D] sigma_between;    // between-treatment SDs
  cholesky_factor_corr[D] Lcorr_between;  // between-treatment correlation
}

transformed parameters {
  matrix[D, D] L_Sigma;                   // Cholesky of within-subject covariance
  matrix[D, D] L_Omega;                   // Cholesky of between-treatment covariance
  array[K] vector[D] mu;                  // treatment-specific mean vectors

  L_Sigma = diag_pre_multiply(sigma_within, Lcorr_within);
  L_Omega = diag_pre_multiply(sigma_between, Lcorr_between);

  for (k in 1:K) {
    mu[k] = mu0 + L_Omega * z_mu[k];
  }
}

model {
  // hyperpriors
  mu0 ~ normal(0, 5);

  sigma_within ~ normal(0, 2);
  sigma_between ~ normal(0, 2);

  Lcorr_within ~ lkj_corr_cholesky(2);
  Lcorr_between ~ lkj_corr_cholesky(2);

  // non-centered treatment effects
  for (k in 1:K) {
    z_mu[k] ~ normal(0, 1);
  }

  // likelihood
  for (i in 1:N) {
    to_vector(y[i]) ~ multi_normal_cholesky(mu[trt[i]], L_Sigma);
  }
}

generated quantities {
  cov_matrix[D] Sigma;
  cov_matrix[D] Omega;
 
  Sigma = multiply_lower_tri_self_transpose(L_Sigma);
  Omega = multiply_lower_tri_self_transpose(L_Omega);
 
}


