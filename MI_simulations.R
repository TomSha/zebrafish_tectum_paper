plot_mean_v_sd <- function(fits){

    mu <- fits[1,,]
    sig <- fits[2,,]
    n_stim <- ncol(mu)

    cols <- c(rgb(1, 0, 0, 0.3), rgb(0, 1, 0, 0.3), rgb(0, 0, 1, 0.3), rgb(1, 1, 0, 0.3))

    plot(log(sig) ~ log(mu), 
            col = "white",
            xlab = "log(mean)",
            ylab = "log(sd)")

    for (i in 1 : n_stim){
        points(log(sig[,i]) ~ log(mu[,i]), 
            pch = 19, 
            col = cols[i])
    }
}

estimate_gamma <- function(fits){

    mu <- fits[1,,]
    n_stim <- ncol(mu)
    shape_params <- rep(NA, n_stim)
    rate_params <- rep(NA, n_stim)

    for (i in 1 : n_stim){
        params = fitdistr(mu[mu[,i] > 0, i], "gamma")
        shape_params[i] <- params[[1]][1]
        rate_params[i] <- params[[1]][2]
    }

    return(list("shape" = shape_params, "rate" = rate_params))

}

estimate_lm <- function(fits){

    mu <- fits[1,,]
    sig <- fits[2,,]
    n_stim <- ncol(mu)
    sigma_params <- rep(NA, n_stim)
    intercept_params <- rep(NA, n_stim)
    coefficient_params <- rep(NA, n_stim)

    for (i in 1 : n_stim){
        index <- (sig[, i] > 0) & (mu[, i] > 0)
        params <- lm(log(sig[index, i]) ~ log(mu[index, i]))
        
        sigma_params[i] <- summary(params)$sigma
        intercept_params[i] <- summary(params)$coefficients[1, 1]
        coefficient_params[i] <- summary(params)$coefficients[2, 1]
    }

    return(list("sigma" = sigma_params, 
                "intercept" = intercept_params, 
                "coefficient" = coefficient_params)
                )

}

simulate_MI <- function(gamma_params, lm_params, n_cells = 2463, n_stim = 4, resample = F, emp_means = NA){

    n_stim <- length(gamma_params$shape) 
    mu_mat <- matrix(nrow = n_cells, ncol = n_stim)
    sd_mat <- matrix(nrow = n_cells, ncol = n_stim)
    MI_vec <- rep(NA, n_cells)
    signal_vec <- rep(NA, n_cells)

    for (i in 1 : n_cells){

        if(resample){
            
            emp_means_pos <- emp_means
            emp_means_pos[emp_means_pos < 0] <- 0.0001
            mu_mat[i, ] <- emp_means_pos[sample(1 : nrow(emp_means_pos), size = 1),]

        }else{
            mu_mat[i, ] <- rgamma(n = n_stim, 
                                shape = gamma_params$shape, 
                                rate = gamma_params$rate)

        }

        sd_mat[i, ] <- exp(rnorm(n = n_stim, 
                                 mean = lm_params$intercept + log(mu_mat[i, ]) * lm_params$coefficient, 
                                 sd = lm_params$sigma))

        prob <- cbind(mu_mat[i,], sd_mat[i, ])
        MI_vec[i] <- MI(prob)

        signal_vec[i] <- max(mu_mat[i, ])
    }

    cor_val <- cor(MI_vec, signal_vec)

    return(list("MI" = MI_vec, "signal" = signal_vec, "cor" = cor_val))


}


source("constants.R")
source("functions/3_MI_functions.R")

prefix <- prefix_list[1]
data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep = "")

MI_grat <- read.table(paste(data_directory, "MI_grat.dat", sep = ""))[,1]
MI_dot <- read.table(paste(data_directory, "MI_dot.dat", sep = ""))[,1]

fits <- readRDS(paste(data_directory, "response_parametric_fits.RDS", sep=""))
dot_fits <- fits$dot
grat_fits <- fits$grat

plot_mean_v_sd(dot_fits)
gamma_params <- estimate_gamma(dot_fits)
lm_params <- estimate_lm(dot_fits)
MI_sims <- simulate_MI(gamma_params, lm_params, resample = T, emp_means = dot_fits[1,,])

MI_cor_ind <- rep(NA, 50)

for (i in 1 : 50){
    cat(i,"\r")
    MI_cor_ind[i] <- simulate_MI(gamma_params, lm_params)$cor
}


MI_cor_resample <- rep(NA, 50)

for (i in 1 : 50){
    cat(i,"\r")
    MI_cor_resample[i] <- simulate_MI(gamma_params, lm_params, resample = T, emp_means = dot_fits[1,,])$cor
}