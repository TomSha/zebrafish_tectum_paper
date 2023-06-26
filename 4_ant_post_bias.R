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
	ldot_model_n_list <- vector("list", n_exps)

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
#		write.table(xy_bins_list[[i]], paste(data_directory, "xy_bins.dat", sep=""), col.names=F, row.names=F)

	}


	# Find which neurons have non_zero MI (MI_stim_thresh_list) and are visually responsive (NCC_thresh_list)
	dot_thresh_list <- mapply(function(MI_thresh, NCC_thresh) MI_thresh & NCC_thresh, MI_dot_thresh_list, NCC_thresh_list)
	grat_thresh_list <- mapply(function(MI_thresh, NCC_thresh) MI_thresh & NCC_thresh, MI_grat_thresh_list, NCC_thresh_list)

	
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

	model_output_dot <- run_multi_model(dat_dot)
	model_output_grat <- run_multi_model(dat_grat)

	dat <- list(dot = dat_dot, grat = dat_grat)
	model_output <- list(dot = model_output_dot, grat = model_output_grat)

	if(!dir.exists(paste(main_directory,"info_analysis/",sep=""))){
		dir.create(paste(main_directory,"info_analysis/",sep=""))
	}
#	saveRDS(model_output,paste(main_directory,"info_analysis/",model_name,".RDS",sep=""))
#	saveRDS(dat,paste(main_directory,"info_analysis/",model_name,"_dat.RDS",sep=""))
}

















#model name = multi_model_subtype (monocular), multi_model_marg_EF (binocular, eyes free), multi_model_marg_EA (binocular, eyes agarose)
save_multi_model_subtypes<-function(model_name,NCC_thresh=T){

	if(NCC_thresh){
		MI_dot_thresh_list<-mapply(function(MI_thresh,NCC_thresh) MI_thresh & NCC_thresh, MI_dot_thresh_list,NCC_thresh_list)
		MI_grat_thresh_list<-mapply(function(MI_thresh,NCC_thresh) MI_thresh & NCC_thresh, MI_grat_thresh_list,NCC_thresh_list)
		MI_ldot_thresh_list<-mapply(function(MI_thresh,NCC_thresh) MI_thresh & NCC_thresh, MI_ldot_thresh_list,NCC_thresh_list)
	}

	#Find distriubtion of neuronal subtypes by bin for dot stimuli
	max_subtype <- max(sapply(dot_model_n_list,max))

	dot_subtype_thresh<-mapply(function(subtype,MI_thresh) sapply(1:max_subtype, function(x) subtype==x & MI_thresh),dot_model_n_list,MI_dot_thresh_list)  
	bins_dot_count<-mapply(function(subtype,xy_bins,NCC_thresh) apply(subtype,2,function(subtype) extract_bins(xy_bins,subtype,NCC_thresh)),dot_subtype_thresh,xy_bins_list,NCC_thresh_list,SIMPLIFY=F)
	bins_dot_count<-simplify2array(bins_dot_count)
	bins_dot_count<-aperm(bins_dot_count,c(1,3,2))

	
	sample_size<-apply(bins_dot_count,c(2,3),sum)
	n_exp<-dim(bins_dot_count)[2]
	n_bin<-dim(bins_dot_count)[1]
	n_marg<-dim(bins_dot_count)[3]
	dat_dot<-list(bins_count=bins_dot_count,sample_size=sample_size,n_exp=n_exp,n_bin=n_bin,n_marg=n_marg)


	#Find distriubtion of neuronal subtypes by bin for grat stimuli
	max_subtype <- max(sapply(grat_model_n_list,max))

	grat_subtype_thresh<-mapply(function(subtype,MI_thresh) sapply(1:max_subtype, function(x) subtype==x & MI_thresh),grat_model_n_list,MI_grat_thresh_list)  
	bins_grat_count<-mapply(function(subtype,xy_bins,NCC_thresh) apply(subtype,2,function(subtype) extract_bins(xy_bins,subtype,NCC_thresh)),grat_subtype_thresh,xy_bins_list,NCC_thresh_list,SIMPLIFY=F)
	bins_grat_count<-simplify2array(bins_grat_count)
	bins_grat_count<-aperm(bins_grat_count,c(1,3,2))

	sample_size<-apply(bins_grat_count,c(2,3),sum)
	n_exp<-dim(bins_grat_count)[2]
	n_bin<-dim(bins_grat_count)[1]
	n_marg<-dim(bins_grat_count)[3]
	dat_grat<-list(bins_count=bins_grat_count,sample_size=sample_size,n_exp=n_exp,n_bin=n_bin,n_marg=n_marg)


	#Find distriubtion of neuronal subtypes by bin for ldot stimuli
	max_subtype <- max(sapply(ldot_model_n_list,max))

	ldot_subtype_thresh<-mapply(function(subtype,MI_thresh) sapply(1:max_subtype, function(x) subtype==x & MI_thresh),ldot_model_n_list,MI_ldot_thresh_list)  
	bins_ldot_count<-mapply(function(subtype,xy_bins,NCC_thresh) apply(subtype,2,function(subtype) extract_bins(xy_bins,subtype,NCC_thresh)),ldot_subtype_thresh,xy_bins_list,NCC_thresh_list,SIMPLIFY=F)
	bins_ldot_count<-simplify2array(bins_ldot_count)
	bins_ldot_count<-aperm(bins_ldot_count,c(1,3,2))

	
	sample_size<-apply(bins_ldot_count,c(2,3),sum)
	n_exp<-dim(bins_ldot_count)[2]
	n_bin<-dim(bins_ldot_count)[1]
	n_marg<-dim(bins_ldot_count)[3]
	dat_ldot<-list(bins_count=bins_ldot_count,sample_size=sample_size,n_exp=n_exp,n_bin=n_bin,n_marg=n_marg)


	# run model on dot and grating data
	model_output_dot<-run_multi_model(dat_dot)
	model_output_grat<-run_multi_model(dat_grat)
	model_output_ldot<-run_multi_model(dat_ldot)
	model_output<-list(dot=model_output_dot,grat=model_output_grat,ldot=model_output_ldot)
	#model_output<-list(grat=model_output_grat)

	dat<-list(dot=dat_dot,grat=dat_grat,ldot=dat_ldot)

#	dat<-list(grat=dat_grat)

	saveRDS(model_output,paste(main_directory,"info_analysis/",model_name,".RDS",sep=""))
	saveRDS(dat,paste(main_directory,"info_analysis/",model_name,"_dat.RDS",sep=""))
}



save_hier_lin_reg <- function(model_name){

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

	dat_dot <- list(N = N, Y = MI_dot, N_group = N_group, X = AP_dot_bias, group = groups)
	writeLines(readLines("hier_lin_reg.stan"))
	mod_fit_dot <- stan(file = "hier_lin_reg.stan", data = dat_dot, refresh = 150, iter = 5000, chains = 4)



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

	dat_grat <- list(N = N, Y = MI_grat, N_group = N_group, X = AP_grat_bias, group = groups)
	writeLines(readLines("hier_lin_reg.stan"))
	mod_fit_grat <- stan(file = "hier_lin_reg.stan", data = dat_grat, refresh = 150, iter = 5000, chains = 4)

	dat <- list(dot = dat_dot, grat = dat_grat)
	model_output <- list(dot = mod_fit_dot, grat = mod_fit_grat)

	saveRDS(model_output, paste(main_directory, "info_analysis/", model_name, ".RDS", sep = ""))
	saveRDS(dat, paste(main_directory, "info_analysis/", model_name, "_dat.RDS", sep = ""))
}


