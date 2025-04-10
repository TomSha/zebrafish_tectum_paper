data {
  int<lower=1> N_obs;
  int<lower=1> Dim;
  matrix[N_obs, Dim] x_obs;
  vector[N_obs] y_obs;
}

parameters {
  real<lower=0> rho1;
  real<lower=0> rho2;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

model {
  matrix[N_obs, N_obs] cov;
  matrix[N_obs, N_obs] L_cov;
  
  for(i in 1:N_obs) {
   cov[i,i] = square(sigma) + square(alpha);
    
    for(j in (i+1):N_obs) {
     cov[i,j] = square(alpha) * exp(-1/(2*square(rho1)) * square(x_obs[i,1] - x_obs[j,1]) -1/(2*square(rho2)) * square(x_obs[i,2] - x_obs[j,2]));	
     cov[j,i] = cov[i,j];
    }
  } 
   
                             
  L_cov  = cholesky_decompose(cov);

  rho1 ~ inv_gamma(10, 1000); 
  rho2 ~ inv_gamma(10, 1000); 
  alpha ~ exponential(1);
  sigma ~ exponential(1);

  y_obs ~ multi_normal_cholesky(rep_vector(0, N_obs), L_cov);
}
