library(bayestestR)

source("functions/4_ant_post_bias_functions.R")
source("models/hier_multinomial.R")


AP_bias <- function(prefix_list){

	n_exps <- length(prefix_list)

	xy_bins_list <- vector("list", n_exps)

	MI_dot_list <- vector("list", n_exps)
	MI_grat_list <- vector("list", n_exps)

	MI_dot_thresh_list <- vector("list", n_exps)
	MI_grat_thresh_list <- vector("list", n_exps)

	dot_model_n_list <- vector("list", n_exps)
	grat_model_n_list <- vector("list", n_exps)

	NCC_thresh_list <- vector("list", n_exps)

	for(i in 1:n_exps){
		prefix <- prefix_list[i]
		data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep="")
		xy_r  <-  read.table(paste(main_directory,  prefix,  "/xy_warped_2.dat",  sep=""))[, 1:2]
	
		MI_dot_list[[i]] <- read.table(paste(data_directory, "MI_dot.dat", sep=""))[, 1]
		MI_grat_list[[i]] <- read.table(paste(data_directory, "MI_grat.dat", sep=""))[, 1]
		
		MI_grat_thresh_list[[i]] <- read.table(paste(data_directory, "MI_grat_thresh.dat", sep=""))[, 1]
		MI_dot_thresh_list[[i]] <- read.table(paste(data_directory, "MI_dot_thresh.dat", sep=""))[, 1]

		NCC_thresh_list[[i]] <- read.table(paste(data_directory, "NCC_thresh.dat", sep=""))[, 1]

		# Here we split the neurons up into 10 equally spaces bins along the AP axis
		xy_bins_list[[i]] <- calculate_bins(xy_r)
		write.table(xy_bins_list[[i]], paste(data_directory, "xy_bins.dat", sep=""), col.names=F, row.names=F)

	}

	
	# dot data
	# get normalised number of neurons per bin
	dot_count <- mapply(norm_cell_num, xy_bins_list, MI_dot_thresh_list, NCC_thresh_list)
	dot_count <- array(dot_count, dim = c(dim(dot_count), 1))
	# neuron number per experiment
	sample_size <- apply(dot_count, c(2, 3), sum)
	# number of experiments 
	n_exp <- dim(dot_count)[2]
	# number of bins
	n_bin <- dim(dot_count)[1]
	# number of subtypes
	n_marg <- dim(dot_count)[3]

	dat_dot <- list(bins_count = dot_count, sample_size = sample_size, n_exp = n_exp, n_bin = n_bin, n_marg = n_marg)

	# grat data
	# get normalised number of neurons per bin
	grat_count <- mapply(norm_cell_num, xy_bins_list, MI_grat_thresh_list, NCC_thresh_list)
	grat_count <- array(grat_count, dim = c(dim(grat_count), 1))
	# neuron number per experiment
	sample_size <- apply(grat_count, c(2, 3), sum)
	# number of experiments 
	n_exp <- dim(grat_count)[2]
	# number of bins
	n_bin <- dim(grat_count)[1]
	# number of subtypes
	n_marg <- dim(grat_count)[3]

	dat_grat <- list(bins_count = grat_count, sample_size = sample_size, n_exp = n_exp, n_bin = n_bin, n_marg = n_marg)

	# run the hierarchical multinomial models
	model_output_dot <- run_multi_model(dat_dot)
	model_output_grat <- run_multi_model(dat_grat)

	dat <- list(dot = dat_dot, grat = dat_grat)
	model_output <- list(dot = model_output_dot, grat = model_output_grat)

	# calculate AP bias index
	AP_bias <- calculate_AP_bias(model_output)

	if(!dir.exists(paste(main_directory,"info_analysis/",sep=""))){
		dir.create(paste(main_directory,"info_analysis/",sep=""))
	}

#	saveRDS(model_output,paste(main_directory,"info_analysis/multi_model.RDS",sep=""))
	saveRDS(dat, paste(main_directory,"info_analysis/multi_model_dat.RDS",sep=""))
	saveRDS(AP_bias, paste(main_directory,"info_analysis/ant_post_bias_all.RDS",sep=""))
}

