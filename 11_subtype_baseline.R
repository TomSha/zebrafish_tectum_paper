library(rstan)

save_hier_lin_reg <- function(prefix_list, model_name = "subtype_snr"){

	# Load data

	xy_bins_list<-vector("list",length(prefix_list))

	MI_dot_list<-vector("list",length(prefix_list))
	MI_grat_list<-vector("list",length(prefix_list))

	MI_dot_thresh_list<-vector("list",length(prefix_list))
	MI_grat_thresh_list<-vector("list",length(prefix_list))
	MI_both_thresh_list<-vector("list",length(prefix_list))

	dot_model_n_list<-vector("list",length(prefix_list))
	grat_model_n_list<-vector("list",length(prefix_list))

	NCC_thresh_list<-vector("list",length(prefix_list))

    snr_list<-vector("list", length(prefix_list))

	for(i in 1:length(prefix_list)){
		prefix<-prefix_list[i]
		data_directory<-paste(main_directory,prefix,"/suite2p/combined/",sep="")
		xy_bins_list[[i]]<-read.table(paste(data_directory,"xy_bins.dat",sep=""))[,1]
		
		MI_dot_list[[i]]<-read.table(paste(data_directory,"MI_dot.dat",sep=""))[,1]
		MI_grat_list[[i]]<-read.table(paste(data_directory,"MI_grat.dat",sep=""))[,1]
		
		MI_grat_thresh_list[[i]]<-read.table(paste(data_directory,"MI_grat_thresh.dat",sep=""))[,1]
		MI_dot_thresh_list[[i]]<-read.table(paste(data_directory,"MI_dot_thresh.dat",sep=""))[,1]
		MI_both_thresh_list[[i]]<-read.table(paste(data_directory,"MI_both_thresh.dat",sep=""))[,1]

		dot_model_n_list[[i]]<-read.table(paste(data_directory,"dot_model_number.dat",sep=""))$V1
		grat_model_n_list[[i]]<-read.table(paste(data_directory,"grat_model_number.dat",sep=""))$V1

		NCC_thresh_list[[i]]<-read.table(paste(data_directory,"NCC_thresh.dat",sep=""))[,1]

        baseline_noise<-read.table(paste(data_directory,"baseline_noise.dat",sep=""))[,1]
        baseline_fluor<-read.table(paste(data_directory,"baseline_fluor.dat",sep=""))[,1]
        snr_list[[i]]<-baseline_fluor / baseline_noise

	}


	# Get the average MI per subtype 
	MI_dot_thresh_list <- mapply(function(x, y) x | y, MI_dot_thresh_list, MI_both_thresh_list)
	MI_grat_thresh_list <- mapply(function(x, y) x | y, MI_grat_thresh_list, MI_both_thresh_list)

	MI_dot <- lapply(1:23, function(x) mapply(function(MI,thresh,NCC,mod) MI[thresh & NCC & mod == x], MI_dot_list, MI_dot_thresh_list, NCC_thresh_list, dot_model_n_list))
	MI_dot <- sapply(MI_dot, function(x) sapply(x, mean, na.rm = T))
	MI_dot <- t(MI_dot[,inc])
	MI_dot <- c(MI_dot)

	snr_dot <- lapply(1:23, function(x) mapply(function(snr,thresh,NCC,mod) snr[thresh & NCC & mod == x], snr_list, MI_dot_thresh_list, NCC_thresh_list, dot_model_n_list))
	snr_dot <- sapply(snr_dot, function(x) sapply(x, mean, na.rm = T))
	snr_dot <- t(snr_dot[,inc])
	snr_dot <- c(snr_dot)
    


	N_group <- length(prefix_list)
	N_cl <- length(inc)
	groups <- rep(1:N_group, each = N_cl)

	thresh <- !is.na(MI_dot) & !is.na(snr_dot)

	MI_dot <- MI_dot[thresh]
	snr_dot <- snr_dot[thresh]
	groups <- groups[thresh]

	N <- length(MI_dot)
	snr_mean <- mean(snr_dot)
	snr_norm <- snr_dot - snr_mean

	# run hierarchical model on dot data
	dat_dot <- list(N = N, Y = MI_dot, N_group = N_group, X = snr_norm, group = groups, X_mean = snr_mean)
	writeLines(readLines("models/hier_lin_reg.stan"))
	mod_fit_dot <- stan(file = "models/hier_lin_reg.stan", data = dat_dot, refresh = 150, iter = 5000, chains = 4)


	# Get the average MI per subtype 

	MI_grat<-lapply(1:23, function(x) mapply(function(MI,thresh,NCC,mod) MI[thresh & NCC & mod == x], MI_grat_list, MI_grat_thresh_list, NCC_thresh_list, grat_model_n_list))
	MI_grat <- sapply(MI_grat, function(x) sapply(x, mean, na.rm = T))
	MI_grat <- t(MI_grat[,inc])
	MI_grat <- c(MI_grat)

	snr_grat <- lapply(1:23, function(x) mapply(function(snr,thresh,NCC,mod) snr[thresh & NCC & mod == x], snr_list, MI_grat_thresh_list, NCC_thresh_list, grat_model_n_list))
	snr_grat <- sapply(snr_grat, function(x) sapply(x, mean, na.rm = T))
	snr_grat <- t(snr_grat[,inc])
	snr_grat <- c(snr_grat)


	N_group <- length(prefix_list)
	N_cl <- length(inc)
	groups <- rep(1:N_group, each = N_cl)

	thresh <- !is.na(MI_grat) & !is.na(snr_grat)

	MI_grat <- MI_grat[thresh]
    snr_grat <- snr_grat[thresh]
	groups <- groups[thresh]

	N <- length(MI_grat)
	snr_mean <- mean(snr_grat)
	snr_norm <- snr_grat - snr_mean

	# run hierarchical model on grat data
	dat_grat <- list(N = N, Y = MI_grat, N_group = N_group, X = snr_norm, group = groups, X_mean = snr_mean)
	writeLines(readLines("models/hier_lin_reg.stan"))
	mod_fit_grat <- stan(file = "models/hier_lin_reg.stan", data = dat_grat, refresh = 150, iter = 5000, chains = 4)

	dat <- list(dot = dat_dot, grat = dat_grat)
	model_output <- list(dot = mod_fit_dot, grat = mod_fit_grat)

	saveRDS(model_output, paste(main_directory, "info_analysis/", model_name, ".RDS", sep = ""))
	saveRDS(dat, paste(main_directory, "info_analysis/", model_name, "_dat.RDS", sep = ""))
}


