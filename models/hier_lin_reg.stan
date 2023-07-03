data {
	int<lower=1> N;
	vector[N] Y;
	int<lower = 1> N_group;
	vector[N] X;
	int<lower = 1, upper = N_group> group[N];
	}
	
parameters {
	real<lower = 0> sigma;
	vector<lower  = 0>[2] sigma_group;
	real B0;
	real B1;
	matrix[2, N_group] z_B;
	cholesky_factor_corr[2] L_B;
	}

transformed parameters {
	matrix[N_group, 2] B;
	B = (diag_pre_multiply(sigma_group, L_B) * z_B)';
}

model {
	B0 ~ normal(1,5);
	B1 ~ normal(0,5);
	sigma ~ exponential(2);
	sigma_group ~ exponential(2);
	L_B ~ lkj_corr_cholesky(1);
	to_vector(z_B) ~ std_normal();

	Y ~ normal(B0 + B[group, 1] + X .* (B1 + B[group, 2]), sigma);
}


generated quantities {

	corr_matrix[2] rho_b = L_B * L_B';
	vector[N_group] effect_by_group;

	for(i in 1:N_group){
		effect_by_group[i] = B1 + B[i, 2];
	}

}



	

