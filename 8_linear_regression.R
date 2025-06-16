library(rstan)

save_hier_lin_reg <- function(prefix_list){

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

	}

	AP_bias <- readRDS(paste(main_directory, "info_analysis/ant_post_bias_subtype.RDS", sep = ""))


	# Get the average MI per subtype 
	MI_dot_thresh_list <- mapply(function(x, y) x | y, MI_dot_thresh_list, MI_both_thresh_list)
	MI_grat_thresh_list <- mapply(function(x, y) x | y, MI_grat_thresh_list, MI_both_thresh_list)

	MI_dot <- lapply(1:23, function(x) mapply(function(MI,thresh,NCC,mod) MI[thresh & NCC & mod == x], MI_dot_list, MI_dot_thresh_list, NCC_thresh_list, dot_model_n_list))
	MI_dot <- sapply(MI_dot, function(x) sapply(x, mean))
	MI_dot <- t(MI_dot[,inc])
	MI_dot <- c(MI_dot)

	AP_dot_bias <- t(AP_bias$dot[,inc])
	AP_dot_bias <- c(AP_dot_bias)

	N_group <- length(prefix_list)
	N_cl <- length(inc)
	groups <- rep(1:N_group, each = N_cl)

	thresh <- !is.na(MI_dot)

	MI_dot <- MI_dot[thresh]
	AP_dot_bias <- AP_dot_bias[thresh]
	groups <- groups[thresh]

	N <- length(AP_dot_bias)

	# run hierarchical model on dot data
	dat_dot <- list(N = N, Y = MI_dot, N_group = N_group, X = AP_dot_bias, group = groups)
	writeLines(readLines("models/hier_lin_reg.stan"))
	mod_fit_dot <- stan(file = "models/hier_lin_reg.stan", data = dat_dot, refresh = 150, iter = 5000, chains = 4)


	# Get the average MI per subtype 

	MI_grat<-lapply(1:23, function(x) mapply(function(MI,thresh,NCC,mod) MI[thresh & NCC & mod == x], MI_grat_list, MI_grat_thresh_list, NCC_thresh_list, grat_model_n_list))
	MI_grat <- sapply(MI_grat, function(x) sapply(x, mean))
	MI_grat <- t(MI_grat[,inc])
	MI_grat <- c(MI_grat)

	AP_grat_bias <- t(AP_bias$grat[,inc])
	AP_grat_bias <- c(AP_grat_bias)

	N_group <- length(prefix_list)
	N_cl <- length(inc)
	groups <- rep(1:N_group, each = N_cl)

	thresh <- !is.na(MI_grat)

	MI_grat <- MI_grat[thresh]
	AP_grat_bias <- AP_grat_bias[thresh]
	groups <- groups[thresh]

	N <- length(AP_grat_bias)

	# run hierarchical model on grat data
	dat_grat <- list(N = N, Y = MI_grat, N_group = N_group, X = AP_grat_bias, group = groups)
	writeLines(readLines("models/hier_lin_reg.stan"))
	mod_fit_grat <- stan(file = "models/hier_lin_reg.stan", data = dat_grat, refresh = 150, iter = 5000, chains = 4)

	dat <- list(dot = dat_dot, grat = dat_grat)
	model_output <- list(dot = mod_fit_dot, grat = mod_fit_grat)

	saveRDS(model_output, paste(main_directory, "info_analysis/", model_name, ".RDS", sep = ""))
	saveRDS(dat, paste(main_directory, "info_analysis/", model_name, "_dat.RDS", sep = ""))
}


