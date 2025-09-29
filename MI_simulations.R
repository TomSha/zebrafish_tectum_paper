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




#####################



source("constants.R")
prefix <- prefix_list[1]
data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep = "")
bl_noise <- read.table(paste(data_directory, "baseline_noise.dat", sep = ""))$V1
bl_f <- read.table(paste(data_directory, "baseline_fluor.dat", sep =""))$V1
MI_dot <- read.table(paste(data_directory, "MI_grat.dat", sep = ""))[,1]
norm_fits <- readRDS(paste(data_directory, "response_parametric_fits.RDS", sep=""))
xy <- read.table(paste(data_directory, "xy.dat", sep = ""))
z <- xy[,3]

n_cells <- length(bl_noise)
#fits <- norm_fits$dot_fits
fits <- norm_fits$grat_fits
signal <- rep(NA, n_cells)

for (i in 1 : n_cells){
    signal[i] <- max(fits[1, i, ])
}

snr_old <- signal / bl_noise
snr_new <- bl_f / bl_noise

X11();plot(log(MI_dot) ~ log(snr_new))

cor(MI_dot, snr_new, use = "complete")
plot(MI_dot ~ signal)
plot(MI_dot ~ bl_noise)

bl_norm <- (bl_noise - mean(bl_noise, na.rm = T)) / sd(bl_noise, na.rm = T)
ind <- bl_norm < 3
plot(MI_dot[ind] ~ bl_noise[ind])
abline(lm(MI_dot[ind] ~ bl_noise[ind]))

plot(signal ~ bl_noise)

mod <- lm(MI_dot ~ snr)
mod2 <- lm(signal ~ bl_noise)
summary(lm(MI_dot ~ signal + bl_noise))
plot(lm(MI_dot ~ bl_noise))

response_var <- fits[2,,]
response_mean <- fits[1,,]

X11();plot(log(response_var) ~ log(response_mean))
a = 4
#plot(log(response_var[,a]) ~ log(response_mean[,a]))
lm(log(response_var[,a]) ~ log(response_mean[,a]))
mod3 <- lm(log(c(response_var)) ~ log(c(response_mean)))
plot(mod3)




alpha = rep(NA, n_cells)

for (i in 1:n_cells){
    alpha[i] = lm(log(response_var[i,]) ~ log(response_mean[i,]))$coefficients[2]
}


n_stim <- 4

noise <- rexp(n = n_cells)
mu <- rexp(n_cells * n_stim)
sig <- rexp(n_cells * n_stim, rate = 3)

mu <- matrix(mu, nrow = n_cells, ncol = n_stim) + noise

# two separate things
# 1. Are good recording conditions causing MI? - look at noise/baseline
# 2. Are we missing some neurons that have significant mutual information but can#t detect due to noise? SNR?
# because increase in signal can increase SNR and neurons with higher signal are more likely to have higher MI (possibly show with simulation)
# then correlation might be expectee
# 
index <- c(fits[2,,]) > 0 & c(fits[1,,]) > 0
means <- c(fits[1,,])[index]
stds <- c(fits[2,,])[index]

means_params = fitdistr(means,"gamma")
hist(means, breaks = 100, freq = FALSE)
curve(dgamma(x, shape = means_params$estimate[1], rate = means_params$estimate[2]), add = T, col = "red")

pl_mod <- lm(log(stds) ~ log(means))
pl_beta <- pl_mod$coefficients
pl_sig <- summary(pl_mod)$sigma
plot(log(stds) ~ log(means))
abline(pl_mod, col = "red")


#stds_params = fitdistr(stds, "gamma")
#hist(stds, breaks = 100, freq = FALSE)
#curve(dgamma(x, shape = stds_params$estimate[1], rate = stds_params$estimate[2]), add = T, col = "red")

n_stim <- 4
n_reps <- 10

mu_mat <- matrix(nrow = n_cells, ncol = n_stim)
sd_mat <- matrix(nrow = n_cells, ncol = n_stim)

signal_vec <- rep(NA, n_cells)
MI_vec <-rep(NA, n_cells)
cor_vec <- rep(NA, 100)

for (j in 1: 100){
    for (i in 1 : n_cells){

        resp_mat <- matrix(nrow = n_reps, ncol = n_stim)
    #    mu <- rep(NA, n_stim)
    #    sd <- rep(1, n_stim)

        mu_mat[i, ] <- rgamma(n = 4, shape = means_params$estimate[1], rate = means_params$estimate[2])
        sd_mat[i, ] <- exp(rnorm(n = 4, mean = pl_beta[[1]] + log(mu_mat[i, ]) * pl_beta[[2]], sd = pl_sig))
    #    sd_mat[i, ] <- exp(rnorm(n = 4, mean = log(mu_mat[i, ]) * 2, sd = 0.0001))

        prob <- cbind(mu_mat[i,], sd_mat[i, ])
        MI_val <- MI(prob)

        signal_vec[i] <- max(mu_mat[i, ])
        MI_vec[i] <- MI_val
    }

    cor_vec[j] <- cor(signal_vec, MI_vec)
}

plot(MI_vec ~ signal_vec)
plot(signal ~ signal_vec)
plot(log(c(sd_mat)) ~ log(c(mu_mat)), col = rgb(0, 0, 0, 0.1), pch = 19)
lines(log(stds) ~ log(means), col = rgb(1, 0, 0, 0.1), type = "p")
plot(log(bl_noise) ~ log(signal))
lm(log())

hist(means, freq = FALSE, breaks = 100)
hist(c(mu_mat), freq = FALSE, breaks = 100, add = T, col = rgb(0, 0, 1, 0.5))

hist(stds, freq = FALSE, breaks = 100)
hist(c(sd_mat), freq = FALSE, breaks = 100, add = T, col = rgb(0, 0, 1, 0.5))

hist(signal, freq = FALSE, breaks = 100)
hist(signal_vec, freq = FALSE, breaks = 100, add = T, col = rgb(0, 0, 1, 0.5))

mod <- lm(MI_dot ~ bl_f + bl_noise + z)
X11()
plot(mod)
X11()
plot(MI_dot ~ z)