data {
  int<lower=2> R;                 // number of diagnostic score categories
  array[2 * R] int<lower=0> y;    // joint counts:
                                  // y[1:R] = affected scores 1,...,R
                                  // y[(R+1):(2R)] = unaffected scores 1,...,R
  vector<lower=0>[2 * R] alpha;   // Dirichlet prior parameters
}

parameters {
  simplex[2 * R] p;               // joint probabilities P(D, S)
}

model {
  p ~ dirichlet(alpha);
  y ~ multinomial(p);
}

generated quantities {
  vector[R] theta_star;           // P(S = i | D = affected)
  vector[R] phi_star;             // P(S = i | D = unaffected)

  real p_affected;
  real p_unaffected;
  real auc;

  p_affected = sum(p[1:R]);
  p_unaffected = sum(p[(R + 1):(2 * R)]);

  for (i in 1:R) {
    theta_star[i] = p[i] / p_affected;
    phi_star[i] = p[R + i] / p_unaffected;
  }

  auc = 0;

  for (i in 2:R) {
    auc += theta_star[i] * sum(phi_star[1:(i - 1)]);
  }

  for (i in 1:R) {
    auc += 0.5 * theta_star[i] * phi_star[i];
  }
}
