
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
    p_val <- rep(NA, n_stim)
    std_err <- rep(NA, n_stim)

    for (i in 1 : n_stim){
        index <- (sig[, i] > 0) & (mu[, i] > 0)
        params <- lm(log(sig[index, i]) ~ log(mu[index, i]))
        
        sigma_params[i] <- summary(params)$sigma
        intercept_params[i] <- summary(params)$coefficients[1, 1]
        coefficient_params[i] <- summary(params)$coefficients[2, 1]
        p_val[i] <- summary(params)$coefficients[2, 4]
        std_err[i] <- summary(params)$coefficients[2, 2]
    }

    return(list("sigma" = sigma_params, 
                "intercept" = intercept_params, 
                "coefficient" = coefficient_params,
                "p_val" = p_val,
                "std_err" = std_err)
                )

}

simulate_MI <- function(gamma_params, lm_params, n_cells, n_stim = 4, resample = F, emp_means = NA){

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

plot_mean_v_sd <- function(fits, i, col = rgb(0, 0, 0, 0.3), params = NA, title = "", S = 1.5){

    mu <- fits[1,,]
    sig <- fits[2,,]
    n_stim <- ncol(mu)

    plot(log(sig) ~ log(mu), 
            col = "white",
            xlab = "log(mean)",
            ylab = "log(sd)",
            main = title,
            cex.lab = S,
            cex.axis = S,
            cex.main = S,
            font = 2,
            font.lab = 2)

    points(log(sig[,i]) ~ log(mu[,i]), 
            pch = 19, 
            col = col)

    if (!is.na(params[1])){
        abline(params$intercept[i], params$coefficient[i], col = "blue")

      x_vals <- seq(min(log(mu[,i]), na.rm = TRUE), max(log(mu[,i]), na.rm = TRUE), length.out = 100)
      y_fit <- params$intercept[i] + params$coefficient[i] * x_vals
      
      # ± noise_sd bands
      polygon(c(x_vals, rev(x_vals)),
        c(y_fit + params$sigma[i] * 2, rev(y_fit - params$sigma[i] * 2)),
        col = adjustcolor("black", 0.2), border = NA)

      polygon(c(x_vals, rev(x_vals)),
        c(y_fit + params$std_err[i], rev(y_fit - params$std_err[i])),
        col = adjustcolor("blue", 0.5), border = NA)
    }
}

plot_gamma <- function(fits, i, gamma_params, main = "", S = 1.5){

    hist(fits[1,,i], breaks = 100, freq = FALSE, xlab = "mean", main = main, xlim = c(0, 6),
            cex.lab = S, cex.axis = S, font = 2, font.lab = 2)
    x <- seq(from = 0, to = 10, length.out = 100)
    curve(dgamma(x, shape = gamma_params$shape[i], rate = gamma_params$rate[i]), add = T, col = "blue")

}


source("constants.R")
source("functions/3_MI_functions.R")

prefix <- prefix_list[2]

save_MI_sims <- function(prefix){

    print(prefix)

    data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep = "")
    MI_dot <- read.table(paste(data_directory, "MI_dot.dat", sep = ""))[,1]
    NCC_thresh <- read.table(paste(data_directory, "NCC_thresh.dat", sep = ""))[,1]
    fits <- readRDS(paste(data_directory, "response_parametric_fits.RDS", sep=""))
    

    thresh <- MI_dot > 0
    MI_dot <- MI_dot[thresh & NCC_thresh]
    dot_fits <- fits$dot[,thresh & NCC_thresh, ]
    n_cells <- length(MI_dot)

    gamma_params <- estimate_gamma(dot_fits)
    lm_params <- estimate_lm(dot_fits)
    MI_sims_resample <- simulate_MI(gamma_params, lm_params, n_cells = n_cells, resample = T, emp_means = dot_fits[1,,])
    MI_sims_ind <- simulate_MI(gamma_params, lm_params, n_cells = n_cells, resample = F)

    MI_cor_ind <- rep(NA, 500)

    for (i in 1 : 500){
        cat(i,"\r")
        MI_cor_ind[i] <- simulate_MI(gamma_params, lm_params)$cor
    }

    data_directory <- paste(main_directory,prefix,"/suite2p/combined/",sep="")
    write.table(MI_cor_ind, paste(data_directory, "MI_cor_ind.dat", sep = ""), col.names = F, row.names = F)


    MI_cor_resample <- rep(NA, 500)

    for (i in 1 : 500){
        cat(i,"\r")
        MI_cor_resample[i] <- simulate_MI(gamma_params, lm_params, resample = T, emp_means = dot_fits[1,,])$cor
    }
    data_directory <- paste(main_directory,prefix,"/suite2p/combined/",sep="")
    write.table(MI_cor_resample, paste(data_directory, "MI_cor_resample.dat", sep = ""), col.names = F, row.names = F)
}


plot_fig <- function(S = 1.5){

    set.seed(4)
    par(mar = c(5.1, 5.1, 4.1, 2.1))

    data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep = "")
    MI_dot <- read.table(paste(data_directory, "MI_dot.dat", sep = ""))[,1]
    NCC_thresh <- read.table(paste(data_directory, "NCC_thresh.dat", sep = ""))[,1]
    MI_cor_ind <-read.table(paste(data_directory, "MI_cor_ind.dat", sep = ""))[,1]
    MI_cor_resample <-read.table(paste(data_directory, "MI_cor_resample.dat", sep = ""))[,1]
    fits <- readRDS(paste(data_directory, "response_parametric_fits.RDS", sep=""))

	timelog <- read.table(paste(main_directory,prefix,"/time_epoches.log",sep=""))
    epoch_labels <- get_epoch_labels(timelog)
    dot_labels <- epoch_labels[grep("^DOT", epoch_labels)]
    dot_ord <- order(dot_labels)

    thresh <- MI_dot > 0
    MI_dot <- MI_dot[thresh & NCC_thresh]
    dot_fits <- fits$dot[,thresh & NCC_thresh, dot_ord]

    n_cells <- length(MI_dot)
    signal <- rep(NA, n_cells)

    for (i in 1 : n_cells){
        signal[i] <- max(dot_fits[1, i, ])
    }

    stim <- group_names[c(4, 3, 1, 2)]

    gamma_params <- estimate_gamma(dot_fits)
    lm_params <- estimate_lm(dot_fits)
    MI_sims_resample <- simulate_MI(gamma_params, lm_params, n_cells = n_cells, resample = T, emp_means = dot_fits[1,,])
    MI_sims_ind <- simulate_MI(gamma_params, lm_params, n_cells = n_cells, resample = F)

    lmat <- matrix(c(1:12), ncol = 4, byrow = T) 
    layout(lmat)


    for(i in 1:4){
        title = paste(stim[i], "\n", 
                        "slope: ", 
                        round(lm_params$coefficient[i], 2) , 
                        " SE: ", round(lm_params$std_err[i], 2))

        plot_mean_v_sd(dot_fits, i, params = lm_params, S = S, title = title)
        if (i == 1) title(main = "A)", adj = 0, cex.main = 2)
    }
    for(i in 1:4){
        plot_gamma(dot_fits, i, gamma_params = gamma_params, S = S)
        if (i == 1) title(main = "B)", adj = 0, cex.main = 2)
    }

    plot(MI_sims_resample$MI ~ MI_sims_resample$signal, 
            pch = 19, 
            col = rgb(0, 0, 0, 0.3),
            xlab = "maximum mean response",
            ylab = "MI",
            font.lab = 2,
            font = 2,
            font.main = 2,
            cex.lab = S,
            cex.axis = S,
            cex.main = S,
            main = "simulations"
            )
    title(main = "C)", adj = 0, cex.main = 2)

    plot(MI_dot ~ signal,
            pch = 19, 
            col = rgb(0, 0, 0, 0.3),
            xlab = "maximum mean response",
            ylab = "MI",
            font.lab = 2,
            font = 2,
            font.main = 2,
            cex.lab = S,
            cex.axis = S,
            cex.main = S,
            main = "Data"
            )
    title(main = "D)", adj = 0, cex.main = 2)

    hist(MI_cor_ind, main = "", 
            xlab = "correlation coefficient", 
            font = 2, 
            font.lab = 2, 
            cex.lab = S, cex.axis = S, xlim = c(0.5, 0.7))


    correlation <- cor(MI_dot, signal)
    abline(v = correlation, col = "blue", lwd = 2)
    title(main = "E)", adj = 0, cex.main = 2)
}

#####################
save_paper_figs_sup <-function(){
	user<-Sys.info()[names(Sys.info())=="user"]
	save_direc<-paste("/media/",user,"/Samsung_T5/Work/tectal_subdomain_paper/paper_figs/",sep="")

	tiff(paste(save_direc, "fig_A2.tiff", sep = ""), units = "cm", width = 30, height = 20, res = 300)
    plot_fig()
	dev.off()

}


#source("constants.R")
#prefix <- prefix_list[1]
#data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep = "")
#bl_noise <- read.table(paste(data_directory, "baseline_noise.dat", sep = ""))$V1
#bl_norm_noise <- read.table(paste(data_directory, "baseline_norm_noise.dat", sep = ""))$V1
#bl_f <- read.table(paste(data_directory, "baseline_fluor.dat", sep =""))$V1
#MI_dot <- read.table(paste(data_directory, "MI_grat.dat", sep = ""))[,1]
#norm_fits <- readRDS(paste(data_directory, "response_parametric_fits.RDS", sep=""))
#xy <- read.table(paste(data_directory, "xy.dat", sep = ""))
#z <- xy[,3]
#
#n_cells <- length(bl_noise)
##fits <- norm_fits$dot_fits
#fits <- norm_fits$grat_fits
#signal <- rep(NA, n_cells)
#
#for (i in 1 : n_cells){
#    signal[i] <- max(fits[1, i, ])
#}
#
#snr_old <- signal / bl_norm_noise
#snr_new <- bl_f / bl_noise
#
#X11();plot(MI_dot ~ snr_old)
#X11();plot(log(MI_dot) ~ log(snr_new))
#
#cor(MI_dot, snr_new, use = "complete")
#plot(MI_dot ~ signal)
#plot(MI_dot ~ bl_noise)
#
#bl_norm <- (bl_noise - mean(bl_noise, na.rm = T)) / sd(bl_noise, na.rm = T)
#ind <- bl_norm < 3
#plot(MI_dot[ind] ~ bl_noise[ind])
#abline(lm(MI_dot[ind] ~ bl_noise[ind]))
#
#plot(signal ~ bl_noise)
#
#mod <- lm(MI_dot ~ snr)
#mod2 <- lm(signal ~ bl_noise)
#summary(lm(MI_dot ~ signal + bl_noise))
#plot(lm(MI_dot ~ bl_noise))
#
#response_var <- fits[2,,]
#response_mean <- fits[1,,]
#
#X11();plot(log(response_var) ~ log(response_mean))
#a = 4
##plot(log(response_var[,a]) ~ log(response_mean[,a]))
#lm(log(response_var[,a]) ~ log(response_mean[,a]))
#mod3 <- lm(log(c(response_var)) ~ log(c(response_mean)))
#plot(mod3)
#
#
#
#
#alpha = rep(NA, n_cells)
#
#for (i in 1:n_cells){
#    alpha[i] = lm(log(response_var[i,]) ~ log(response_mean[i,]))$coefficients[2]
#}
#
#
#n_stim <- 4
#
#noise <- rexp(n = n_cells)
#mu <- rexp(n_cells * n_stim)
#sig <- rexp(n_cells * n_stim, rate = 3)
#
#mu <- matrix(mu, nrow = n_cells, ncol = n_stim) + noise
#
## two separate things
## 1. Are good recording conditions causing MI? - look at noise/baseline
## 2. Are we missing some neurons that have significant mutual information but can#t detect due to noise? SNR?
## because increase in signal can increase SNR and neurons with higher signal are more likely to have higher MI (possibly show with simulation)
## then correlation might be expectee
## 
#index <- c(fits[2,,]) > 0 & c(fits[1,,]) > 0
#means <- c(fits[1,,])[index]
#stds <- c(fits[2,,])[index]
#
#means_params = fitdistr(means,"gamma")
#hist(means, breaks = 100, freq = FALSE)
#curve(dgamma(x, shape = means_params$estimate[1], rate = means_params$estimate[2]), add = T, col = "red")
#
#pl_mod <- lm(log(stds) ~ log(means))
#pl_beta <- pl_mod$coefficients
#pl_sig <- summary(pl_mod)$sigma
#plot(log(stds) ~ log(means))
#abline(pl_mod, col = "red")
#
#
##stds_params = fitdistr(stds, "gamma")
##hist(stds, breaks = 100, freq = FALSE)
##curve(dgamma(x, shape = stds_params$estimate[1], rate = stds_params$estimate[2]), add = T, col = "red")
#
#n_stim <- 4
#n_reps <- 10
#
#mu_mat <- matrix(nrow = n_cells, ncol = n_stim)
#sd_mat <- matrix(nrow = n_cells, ncol = n_stim)
#
#signal_vec <- rep(NA, n_cells)
#MI_vec <-rep(NA, n_cells)
#cor_vec <- rep(NA, 100)
#
#for (j in 1: 100){
#    for (i in 1 : n_cells){
#
#        resp_mat <- matrix(nrow = n_reps, ncol = n_stim)
#    #    mu <- rep(NA, n_stim)
#    #    sd <- rep(1, n_stim)
#
#        mu_mat[i, ] <- rgamma(n = 4, shape = means_params$estimate[1], rate = means_params$estimate[2])
#        sd_mat[i, ] <- exp(rnorm(n = 4, mean = pl_beta[[1]] + log(mu_mat[i, ]) * pl_beta[[2]], sd = pl_sig))
#    #    sd_mat[i, ] <- exp(rnorm(n = 4, mean = log(mu_mat[i, ]) * 2, sd = 0.0001))
#
#        prob <- cbind(mu_mat[i,], sd_mat[i, ])
#        MI_val <- MI(prob)
#
#        signal_vec[i] <- max(mu_mat[i, ])
#        MI_vec[i] <- MI_val
#    }
#
#    cor_vec[j] <- cor(signal_vec, MI_vec)
#}
#
#plot(MI_vec ~ signal_vec)
#plot(signal ~ signal_vec)
#plot(log(c(sd_mat)) ~ log(c(mu_mat)), col = rgb(0, 0, 0, 0.1), pch = 19)
#lines(log(stds) ~ log(means), col = rgb(1, 0, 0, 0.1), type = "p")
#plot(log(bl_noise) ~ log(signal))
#lm(log())
#
#hist(means, freq = FALSE, breaks = 100)
#hist(c(mu_mat), freq = FALSE, breaks = 100, add = T, col = rgb(0, 0, 1, 0.5))
#
#hist(stds, freq = FALSE, breaks = 100)
#hist(c(sd_mat), freq = FALSE, breaks = 100, add = T, col = rgb(0, 0, 1, 0.5))
#
#hist(signal, freq = FALSE, breaks = 100)
#hist(signal_vec, freq = FALSE, breaks = 100, add = T, col = rgb(0, 0, 1, 0.5))
#
#mod <- lm(MI_dot ~ bl_f + bl_noise + z)
#X11()
#plot(mod)
#X11()
#plot(MI_dot ~ z)
#
#
#
#
#
#
#