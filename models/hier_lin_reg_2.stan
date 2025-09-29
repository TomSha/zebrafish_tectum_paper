data {
  int<lower=1> N;
  vector[N] Y;
  int<lower = 1> N_group;
  vector[N] X1;
  vector[N] X2;
  int<lower = 1, upper = N_group> group[N];
}
	
parameters {
  real<lower = 0> sigma;
  vector<lower = 0>[3] sigma_group;  // intercept + 2 slopes
  real B0;
  real B1;
  real B2;
  matrix[3, N_group] z_B;
  cholesky_factor_corr[3] L_B;
}

transformed parameters {
  matrix[N_group, 3] B;
  B = (diag_pre_multiply(sigma_group, L_B) * z_B)';
}

model {
  // priors
  B0 ~ normal(1, 5);
  B1 ~ normal(0, 5);
  B2 ~ normal(0, 5);
  sigma ~ exponential(2);
  sigma_group ~ exponential(2);
  L_B ~ lkj_corr_cholesky(1);
  to_vector(z_B) ~ std_normal();

  // likelihood
  Y ~ normal(
    B0 + B[group, 1] +
    X1 .* (B1 + B[group, 2]) +
    X2 .* (B2 + B[group, 3]),
    sigma
  );
}

generated quantities {
  corr_matrix[3] rho_b = L_B * L_B';
  vector[N_group] effect1_by_group;
  vector[N_group] effect2_by_group;

  for (i in 1:N_group) {
    effect1_by_group[i] = B1 + B[i, 2];
    effect2_by_group[i] = B2 + B[i, 3];
  }
}
