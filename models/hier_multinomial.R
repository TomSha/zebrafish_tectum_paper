


# changed alpha to unif because the sampler was getting stuck when using a gamma dist. 
# I think if a subtype has no cells in one experiment it doesn't like it
# alpha[i,j] ~ dgamma(0.1,1)

multi_model = "model {

	for (k in 1:n_bin) {
		beta0_base[k] <- 1
	}

	beta0[1:n_bin] ~ ddirch(beta0_base[])

	eta ~ dgamma(2.0, 0.01)   
	
	for(i in 1:n_marg){
		for(j in 1:n_bin){
			alpha[i,j] <- eta * beta0[j]
		}
	}

	for(i in 1:n_exp){
		for(j in 1:n_marg){
			p[1:n_bin,i,j] ~ ddirich(alpha[j,])
			}
		}	

	for(i in 1:n_exp){
			for(j in 1:n_marg){
					bins_count[,i,j] ~ dmulti(p[1:n_bin,i,j],sample_size[i,j])
			}
		}
		 	  
         
}"

run_multi_model <- function(dat){

	varnames = c("p", "alpha")
	burn_in = 1000
	steps = 10000
	thin = 5

	library(rjags)
	fileConn = file("model.tmp")
	writeLines(multi_model, fileConn);
	close(fileConn)

	m = jags.model(file = "model.tmp", data = dat, n.chains = 4);
	update(m,burn_in)
	draw = jags.samples(m, steps, thin = thin, variable.names = varnames)
	return(draw)
}



multi_model_subtype = "model {
#### Global base distribution ####
  for (k in 1:n_bin) {
    beta0_base[k] <- 1
  }
  beta0[1:n_bin] ~ ddirch(beta0_base[])     # global mean across all subtypes

  #### Hyperparameter for how tightly subtype bases cluster around beta0 ####
  eta ~ dgamma(2.0, 0.01)   

  #### Loop over subtypes ####
  for (j in 1:n_marg) {

    # Subtype-specific base vector drawn around global mean
    for (k in 1:n_bin) {
      alpha_beta[j,k] <- eta * beta0[k]
    }
    beta_j[j,1:n_bin] ~ ddirch(alpha_beta[j,])   # subtype-level mean distribution

    # Subtype-specific concentration (controls per-fish variability)
    tau[j] ~ dgamma(2.0, 0.01)

    # Combine concentration and base vector into Dirichlet parameters
    for (k in 1:n_bin) {
      alpha[j,k] <- tau[j] * beta_j[j,k] + 1e-6
    }

    # Mean distribution for subtype j (around beta_j)
    theta_alpha[j,1:n_bin] ~ ddirch(alpha[j,])

    #### Per-experiment multinomial probabilities ####
    for (i in 1:n_exp) {
      p[1:n_bin,i,j] ~ ddirch(theta_alpha[j,])
      bins_count[1:n_bin,i,j] ~ dmulti(p[,i,j], sample_size[i,j])
    }
  }

}"



run_multi_model_subtype <- function(dat){

	varnames = c("beta0_base", "beta0", "beta_j", "tau", "alpha", "theta_alpha", "eta", "p")
	burn_in = 1000
	steps = 10000
	thin = 2

	inits <- function() list(
 	   beta0 = rep(1/dat$n_bin, dat$n_bin),
 	   eta = 10,
 	   tau = rep(10, dat$n_marg)
  )

	library(rjags)
	fileConn = file("model.tmp")
	writeLines(multi_model_subtype, fileConn);
	close(fileConn)

	m = jags.model(file = "model.tmp", data = dat, n.chains = 4, n.adapt = 1000);
	update(m,burn_in)
	draw = jags.samples(m, steps, thin = thin, variable.names = varnames)
	return(draw)
}


#theta_alpha <- apply(model_output_dot$theta_alpha, c(1,2), mean)
#ap_bias = apply(theta_alpha[,1:5],1,sum) - apply(theta_alpha[,6:10],1,sum)
#par(mfrow = c(1,3))
#barplot(ap_bias, ylim = c(-1,1))
#hist(model_output_dot$eta)
#
#theta_alpha <- apply(model_output_grat$theta_alpha, c(1,2), mean)
#ap_bias = apply(theta_alpha[,1:5],1,sum) - apply(theta_alpha[,6:10],1,sum)
#par(mfrow = c(1,3))
#barplot(ap_bias, ylim = c(-1,1))
#hist(model_output_grat$eta)


library(LaplacesDemon)

n_bin <- 10
n_subtype <- 2
n_exp <- 3
beta0_base <- rep(1, n_bin)
beta0 <- rep(NA, n_bin)
alpha_beta <- matrix(nrow = n_subtype, ncol = n_bin)
p <- array(dim = c(n_bin, n_exp, n_subtype))
bins_count <- array(dim = c(n_bin, n_exp, n_subtype))
sample_size <- matrix(nrow = n_exp, ncol = n_subtype, rep(100, n_exp * n_subtype))


beta0[1:n_bin] <- rdirichlet(1, beta0_base)     # global mean across all subtypes
#### Hyperparameter for how tightly subtype bases cluster around beta0 ####
eta <- rgamma(1, 2.0, 0.01)   # mean ≈ 10, controls similarity of beta_j to beta0

#### Loop over subtypes ####
for (j in 1:n_subtype) {
# Subtype-specific base vector drawn around global mean
	for (k in 1:n_bin) {
		alpha_beta[j,k] <- eta * beta0[k]
	}
#### Per-experiment multinomial probabilities ####
	for (i in 1:n_exp) {
		p[1:n_bin,i,j] <- rdirichlet(1, alpha_beta[j,])
		bins_count[1:n_bin,i,j] <- rmultinom(n = 1, prob = p[,i,j], size = sample_size[i,j])
	}
}
