library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

source("functions/4_ant_post_bias_functions.R")
source("functions/5_GP_functions.R")

save_GP_fit <- function(prefix){

	print(prefix)

	# Load data
	data_directory <- paste(main_directory,prefix,"/suite2p/combined/",sep="")
	NCC_thresh <- read.table(paste(data_directory,"NCC_thresh.dat",sep=""))[,1]

	MI_dot_thresh <- read.table(paste(data_directory,"MI_dot_thresh.dat",sep=""))[,1]
	MI_dot <- read.table(paste(data_directory,"MI_dot.dat",sep=""))[,1]

	MI_grat_thresh <- read.table(paste(data_directory,"MI_grat_thresh.dat",sep=""))[,1]
	MI_grat <- read.table(paste(data_directory,"MI_grat.dat",sep=""))[,1]
	
	xy_r <- read.table(paste(main_directory,prefix,"/xy_warped_2.dat",sep=""))

	angle <- 25
	alpha <- angle * pi / 180
	xy_r <- rotation_mat(alpha,xy_r[,1:2])

	dot_pnts_r <- xy_r[MI_dot_thresh & NCC_thresh,]
	grat_pnts_r <- xy_r[MI_grat_thresh & NCC_thresh,]

	MI_dot <- MI_dot[MI_dot_thresh & NCC_thresh]
	MI_grat <- MI_grat[MI_grat_thresh & NCC_thresh]
	
	# Generate X coordinates to make predictions
	marg1_X <- seq(-37, 245, length.out = 100)
	marg2_X <- seq(50, 256, length.out = 100)

	X_predict <- expand.grid(marg1_X,marg2_X)

	# Threshold based on the outline of the tectum
	thresh1 <- X_predict[,2] > X_predict[,1] * -0.7 + 120
	thresh2 <- X_predict[,2] > X_predict[,1] * 0.95 - 55
	thresh3 <- X_predict[,2] < X_predict[,1] * 0.7 + 195
	thresh4 <- X_predict[,2] < X_predict[,1] * -0.5 + 300

	thresh <- thresh1 & thresh2 & thresh3 & thresh4
	X_predict <- X_predict[thresh,]

	writeLines(readLines("fit_GP2d_dat.stan"))
	
	# Run model for dots
	X <- dot_pnts_r
	Y <- MI_dot

	dat <- list(N_obs = nrow(X), Dim = 2, x_obs = X, y_obs = Y)
	fit <- stan(file = "fit_GP2d_dat.stan", data = dat, chains = 4, cores = 4, refresh = 500, iter = 2000)

	samples <- extract(fit)
	n_samps <- length(samples[[1]])

	pred_list <- calc_reg(samples, X, Y, X_predict, n_samps = n_samps)
	
	#saveRDS(fit, paste(data_directory, "GP_stan_dot_fit.RDS", sep = ""))
	#saveRDS(pred_list, paste(data_directory, "GP_stan_dot_predictions.RDS", sep = ""))


	# Run model for gratings
	X <- grat_pnts_r
	Y <- MI_grat

	dat <- list(N_obs = nrow(X), Dim = 2, x_obs = X, y_obs = Y)
	fit <- stan(file = "fit_GP2d_dat.stan", data = dat, chains = 4, cores = 4, refresh = 500, iter = 2000)

	samples <- extract(fit)
	n_samps <- length(samples[[1]])

	pred_list <- calc_reg(samples, X, Y, X_predict, n_samps = n_samps)
	
	#saveRDS(fit, paste(data_directory, "GP_stan_grat_fit.RDS", sep = ""))
	#saveRDS(pred_list, paste(data_directory, "GP_stan_grat_predictions.RDS", sep = ""))

		
}	
