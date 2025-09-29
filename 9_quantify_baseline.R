library(rstan)
source("functions/baseline_functions.R")
source("functions/3_MI_functions.R")


quant_baseline <- function(prefix){
    cat(prefix)

	data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep = "")
    spikes <- read.table(paste(data_directory, "spikes.dat", sep = ""))
	traces <- np$load(paste(data_directory, "traces.npy", sep = ""))

    baseL_index <- calculate_baseline_index(spikes)
#    baseL_index <- as.matrix(read.table(paste(data_directory, "baseL_index.dat", sep = "")))
    bl_quant <- quant_baseline_metrics(traces, baseL_index)


	write.table(baseL_index, paste(data_directory, "baseL_index.dat", sep = ""), row.names = F, col.names = F)
    write.table(bl_quant$noise, paste(data_directory, "baseline_noise.dat", sep = ""), row.names = F, col.names = F)
    write.table(bl_quant$fluor, paste(data_directory, "baseline_fluor.dat", sep = ""), row.names = F, col.names = F)
}


save_snr_hier_lin_reg <- function(prefix_list, model_name = "snr_hier_lm"){

	# Load data

	xy_bins_list<-vector("list",length(prefix_list))

	MI_dot_list<-vector("list",length(prefix_list))
	MI_grat_list<-vector("list",length(prefix_list))

	MI_dot_thresh_list<-vector("list",length(prefix_list))
	MI_grat_thresh_list<-vector("list",length(prefix_list))
	MI_both_thresh_list<-vector("list",length(prefix_list))

	NCC_thresh_list<-vector("list",length(prefix_list))

    snr_list<-vector("list",length(prefix_list))

	for(i in 1:length(prefix_list)){
		prefix<-prefix_list[i]
		data_directory<-paste(main_directory,prefix,"/suite2p/combined/",sep="")
		
		MI_dot_list[[i]]<-read.table(paste(data_directory,"MI_dot.dat",sep=""))[,1]
		MI_grat_list[[i]]<-read.table(paste(data_directory,"MI_grat.dat",sep=""))[,1]
		
		MI_grat_thresh_list[[i]]<-read.table(paste(data_directory,"MI_grat_thresh.dat",sep=""))[,1]
		MI_dot_thresh_list[[i]]<-read.table(paste(data_directory,"MI_dot_thresh.dat",sep=""))[,1]
		MI_both_thresh_list[[i]]<-read.table(paste(data_directory,"MI_both_thresh.dat",sep=""))[,1]

		NCC_thresh_list[[i]]<-read.table(paste(data_directory,"NCC_thresh.dat",sep=""))[,1]
        baseline_noise<-read.table(paste(data_directory,"baseline_noise.dat",sep=""))[,1]
        baseline_fluor<-read.table(paste(data_directory,"baseline_fluor.dat",sep=""))[,1]
        snr_list[[i]]<-baseline_fluor / baseline_noise

	}

    # MI dot vs SNR

	MI_dot <- mapply(function(MI,NCC) MI[NCC], MI_dot_list, NCC_thresh_list)
	snr <- mapply(function(snr,NCC) snr[NCC], snr_list, NCC_thresh_list)

	N_group <- length(prefix_list)
	groups <- sapply(1:N_group, function(x) rep(x, length(MI_dot[[x]])))
    groups <- unlist(groups)


    MI_dot <- unlist(MI_dot)
    snr <- unlist(snr)

    thresh <- !is.na(snr)
    MI_dot <- MI_dot[thresh]
    snr <- snr[thresh]
    groups <- groups[thresh]

	N <- length(MI_dot)
	snr_mean <- mean(snr)
	snr_norm <- snr - snr_mean

	# run hierarchical model on dot data
	dat_dot <- list(N = N, Y = MI_dot, N_group = N_group, X = snr_norm, group = groups, X_mean = snr_mean)
	writeLines(readLines("models/hier_lin_reg.stan"))
	mod_fit_dot <- stan(file = "models/hier_lin_reg.stan", data = dat_dot, refresh = 150, iter = 5000, chains = 4)


    # MI dot (intersection) vs SNR
#	MI_dot_int <- mapply(function(MI,thresh,NCC) MI[thresh & NCC], MI_dot_list, MI_both_thresh_list, NCC_thresh_list)
#	snr <- mapply(function(snr,thresh,NCC) snr[thresh & NCC], snr_list, MI_both_thresh_list, NCC_thresh_list)
#
#	N_group <- length(prefix_list)
#	groups <- sapply(1:N_group, function(x) rep(x, length(MI_dot_int[[x]])))
#    groups <- unlist(groups)
#
#
#    MI_dot_int <- unlist(MI_dot_int)
#    snr <- unlist(snr)
#
#    thresh <- !is.na(snr)
#    MI_dot_int <- MI_dot_int[thresh]
#    snr <- snr[thresh]
#    groups <- groups[thresh]
#
#	N <- length(MI_dot_int)
#	snr_mean <- mean(snr)
#	snr_norm <- snr - snr_mean
#
#	# run hierarchical model on dot data
#	dat_dot_int <- list(N = N, Y = MI_dot_int, N_group = N_group, X = snr_norm, group = groups, X_mean = snr_mean)
#	writeLines(readLines("models/hier_lin_reg.stan"))
#	mod_fit_dot_int <- stan(file = "models/hier_lin_reg.stan", data = dat_dot_int, refresh = 150, iter = 5000, chains = 4)
	#5000 4


    # MI grat vs SNR

	MI_grat <- mapply(function(MI,thresh,NCC) MI[NCC], MI_grat_list, NCC_thresh_list)
	snr <- mapply(function(snr,thresh,NCC) snr[NCC], snr_list, NCC_thresh_list)

	N_group <- length(prefix_list)
	groups <- sapply(1:N_group, function(x) rep(x, length(MI_grat[[x]])))
    groups <- unlist(groups)


    MI_grat <- unlist(MI_grat)
    snr <- unlist(snr)

    thresh <- !is.na(snr)
    MI_grat <- MI_grat[thresh]
    snr <- snr[thresh]
    groups <- groups[thresh]

	N <- length(MI_grat)
	snr_mean <- mean(snr)
	snr_norm <- snr - snr_mean

	# run hierarchical model on grat data
	dat_grat <- list(N = N, Y = MI_grat, N_group = N_group, X = snr_norm, group = groups, X_mean = snr_mean)
	writeLines(readLines("models/hier_lin_reg.stan"))
	mod_fit_grat <- stan(file = "models/hier_lin_reg.stan", data = dat_grat, refresh = 150, iter = 5000, chains = 4)


    # MI grat (intersection) vs SNR
#	MI_grat_int <- mapply(function(MI,thresh,NCC) MI[thresh & NCC], MI_grat_list, MI_both_thresh_list, NCC_thresh_list)
#	snr <- mapply(function(snr,thresh,NCC) snr[thresh & NCC], snr_list, MI_both_thresh_list, NCC_thresh_list)
#
#	N_group <- length(prefix_list)
#	groups <- sapply(1:N_group, function(x) rep(x, length(MI_grat_int[[x]])))
#    groups <- unlist(groups)
#
#
#    MI_grat_int <- unlist(MI_grat_int)
#    snr <- unlist(snr)
#
#    thresh <- !is.na(snr)
#    MI_grat_int <- MI_grat_int[thresh]
#    snr <- snr[thresh]
#    groups <- groups[thresh]
#
#	N <- length(MI_grat_int)
#	snr_mean <- mean(snr)
#	snr_norm <- snr - snr_mean
#
#	# run hierarchical model on grat data
#	dat_grat_int <- list(N = N, Y = MI_grat_int, N_group = N_group, X = snr_norm, group = groups, X_mean = snr_mean)
#	writeLines(readLines("models/hier_lin_reg.stan"))
#	mod_fit_grat_int <- stan(file = "models/hier_lin_reg.stan", data = dat_grat_int, refresh = 150, iter = 5000, chains = 4)


#	dat <- list(dot = dat_dot, dot_int = dat_dot_int, grat = dat_grat, grat_int = dat_grat_int)
#	model_output <- list(dot = mod_fit_dot, dot_int = mod_fit_dot_int, grat = mod_fit_grat, grat_int = mod_fit_grat_int)

	dat <- list(dot = dat_dot, grat = dat_grat)
	model_output <- list(dot = mod_fit_dot, grat = mod_fit_grat)


	saveRDS(model_output, paste(main_directory, "info_analysis/", model_name, ".RDS", sep = ""))
	saveRDS(dat, paste(main_directory, "info_analysis/", model_name, "_dat.RDS", sep = ""))
}