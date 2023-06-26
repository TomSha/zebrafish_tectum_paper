


# changed alpha to unif because the sampler was getting stuck when using a gamma dist. 
# I think if a subtype has no cells in one experiment it doesn't like it
# alpha[i,j] ~ dgamma(0.1,1)

multi_model = "model {
	
	for(i in 1:n_marg){
		for(j in 1:n_bin){
			alpha[i,j] ~ dunif(1,50)
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

