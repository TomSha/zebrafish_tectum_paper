library(rethinking)
calc_r2_bayes <- function(dat, mod_fit){

	params <- extract.samples(mod_fit)
    r2_list <- vector("list", dat$N_group)
    n_samples <- 10000

    for (group in 1 : dat$N_group){
        y <- dat$Y[dat$group == group]
        x <- dat$X[dat$group == group]
        y_hat <- matrix(nrow = n_samples, ncol = length(x))

        for(j in 1:n_samples){
            y_hat[j,] <- params$B0[j] + params$B[j,group,1] + (params$B1[j] + params$B[j,group,2]) * x

        }
        var_fit <- apply(y_hat, 1, var)
        var_resid <- params$sigma^2
        r2_list[[group]] <- var_fit / (var_fit + var_resid)

    }

    return(r2_list)
    


}
