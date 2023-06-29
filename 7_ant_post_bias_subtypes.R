source("functions/4_ant_post_bias_functions.R")

AP_bias_subtypes <- function(prefix_list){
	xy_bins_list <- vector("list",length(prefix_list))

	MI_dot_list <- vector("list",length(prefix_list))
	MI_grat_list <- vector("list",length(prefix_list))

	MI_dot_thresh_list <- vector("list",length(prefix_list))
	MI_grat_thresh_list <- vector("list",length(prefix_list))

	dot_model_n_list <- vector("list",length(prefix_list))
	grat_model_n_list <- vector("list",length(prefix_list))

	NCC_thresh_list <- vector("list",length(prefix_list))

	for(i in 1:length(prefix_list)){
		prefix <- prefix_list[i]
		data_directory <- paste(main_directory,prefix,"/suite2p/combined/",sep="")
		xy_bins_list[[i]] <- read.table(paste(data_directory,"xy_bins.dat",sep=""))[,1]
		
		MI_dot_list[[i]] <- read.table(paste(data_directory,"MI_dot.dat",sep=""))[,1]
		MI_grat_list[[i]] <- read.table(paste(data_directory,"MI_grat.dat",sep=""))[,1]
		
		MI_grat_thresh_list[[i]] <- read.table(paste(data_directory,"MI_grat_thresh.dat",sep=""))[,1]
		MI_dot_thresh_list[[i]] <- read.table(paste(data_directory,"MI_dot_thresh.dat",sep=""))[,1]

		dot_model_n_list[[i]] <- read.table(paste(data_directory,"dot_model_number.dat",sep=""))$V1
		grat_model_n_list[[i]] <- read.table(paste(data_directory,"grat_model_number.dat",sep=""))$V1

		NCC_thresh_list[[i]] <- read.table(paste(data_directory,"NCC_thresh.dat",sep=""))[,1]

	}


	# dot dat
	# get the normalised cell counts per bin for all subtypes
	dot_count <- mapply(norm_cell_num, xy_bins_list, MI_dot_thresh_list, NCC_thresh_list, dot_model_n_list, SIMPLIFY = F)
	dot_count <- simplify2array(dot_count)
	dot_count <- aperm(bins_dot_count, c(1, 3, 2))
	
	sample_size <- apply(bins_dot_count, c(2, 3), sum)
	n_exp <- dim(bins_dot_count)[2]
	n_bin <- dim(bins_dot_count)[1]
	n_marg <- dim(bins_dot_count)[3]
	dat_dot <- list(bins_count = bins_dot_count, sample_size = sample_size, n_exp = n_exp, n_bin = n_bin, n_marg = n_marg)


	# grat dat
	# get the normalised cell counts per bin for all subtypes
	grat_count <- mapply(norm_cell_num, xy_bins_list, MI_grat_thresh_list, NCC_thresh_list, grat_model_n_list, SIMPLIFY = F)
	grat_count <- simplify2array(bins_grat_count)
	grat_count <- aperm(bins_grat_count, c(1, 3, 2))

	sample_size <- apply(bins_grat_count, c(2, 3), sum)
	n_exp <- dim(bins_grat_count)[2]
	n_bin <- dim(bins_grat_count)[1]
	n_marg <- dim(bins_grat_count)[3]
	dat_grat <- list(bins_count = bins_grat_count, sample_size = sample_size, n_exp = n_exp, n_bin = n_bin, n_marg = n_marg)

	# run model on dot and grating data
	model_output_dot <- run_multi_model(dat_dot)
	model_output_grat <- run_multi_model(dat_grat)
	
	model_output <- list(dot = model_output_dot, grat = model_output_grat)

	dat <- list(dot = dat_dot, grat = dat_grat)

	#saveRDS(model_output,paste(main_directory,"info_analysis/multi_model_subtype.RDS",sep=""))
	saveRDS(dat,paste(main_directory,"info_analysis/multi_model_subtype_dat.RDS",sep=""))
}



