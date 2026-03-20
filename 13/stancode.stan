data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y1;
  vector[N] y2;
  vector[N] y3;
}

parameters {
 
  real beta;
  real<lower=0> tau;

  vector[N] eta;

  real nu1;
  real nu2;
  real nu3;

  real lambda2;
  real lambda3;

  real<lower=0> sigma1;
  real<lower=0> sigma2;
  real<lower=0> sigma3;
}

model {
  // Priors
  
  beta  ~ normal(0, 2);
  tau   ~ normal(0, 1);

  nu1 ~ normal(0, 3);
  nu2 ~ normal(0, 3);
  nu3 ~ normal(0, 3);

  lambda2 ~ normal(0, 1.5);
  lambda3 ~ normal(0, 1.5);

  sigma1 ~ normal(0, 1);
  sigma2 ~ normal(0, 1);
  sigma3 ~ normal(0, 1);

  // Structural model
  eta ~ normal(beta * x, tau);

  // Measurement model
  y1 ~ normal(nu1 + eta, sigma1);          // lambda1 fixed to 1
  y2 ~ normal(nu2 + lambda2 * eta, sigma2);
  y3 ~ normal(nu3 + lambda3 * eta, sigma3);
}

// generated quantities {
//   vector[N] y1_rep;
//   vector[N] y2_rep;
//   vector[N] y3_rep;
// 
//   for (i in 1:N) {
//     y1_rep[i] = normal_rng(nu1 + eta[i], sigma1);
//     y2_rep[i] = normal_rng(nu2 + lambda2 * eta[i], sigma2);
//     y3_rep[i] = normal_rng(nu3 + lambda3 * eta[i], sigma3);
//   }
// }

