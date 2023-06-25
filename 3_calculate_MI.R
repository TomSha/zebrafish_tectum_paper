source("functions/3_MI_functions.R")

# Here we calculate the mutual information for the stimulus and response to dots and gratings
calculate_MI <- function(prefix){
	data_directory<-paste(main_directory,prefix,"/suite2p/combined/",sep="")
	timelog<-read.table(paste(main_directory,prefix,"/time_epoches.log",sep=""))
	traces <- np$load(paste(data_directory, "traces_norm.npy", sep = ""))

	responses_split_list <- split_traces(timelog, traces, ac_freq)

	# take the mean response across each epoch
	responses_mean_list <- lapply(responses_split_list, function(x) apply(x, c(1, 3), mean))
	epoch_labels <- get_epoch_labels(timelog)

	# separate out grating and dot responses
	grat_responses_list <- responses_mean_list[grep("GRAT", epoch_labels)]
	dot_responses_list <- responses_mean_list[grep("^DOT", epoch_labels)]

	print("estimating response distributions")

	# estimate P(R|S) for each neuron and stimulus
	grat_fits <- simplify2array(lapply(grat_responses_list, fit_dist))
	dot_fits <- simplify2array(lapply(dot_responses_list, fit_dist))

	print("calculating mutual information")
	
	# calculate the mutual information between the response and the stimuli for dots and gratings
	MI0_grat <- apply(grat_fits, 2, function(x) MI(t(x)))
	MI0_dot <- apply(dot_fits, 2, function(x) MI(t(x)))

	print("calculating shuffled mutual information")

	# shuffle responses and calculate MI to create random noise model	
	MI_grat_shuffle <- calculate_MI_shuffle(grat_responses_list)
	MI_dot_shuffle <- calculate_MI_shuffle(dot_responses_list) 
		
	# subtract the mean random MI from the MI
	MI_grat <- MI0_grat - rowMeans(simplify2array(MI_grat_shuffle))
	MI_dot <- MI0_dot - rowMeans(simplify2array(MI_dot_shuffle))

	# calulate thresholds	
	MI_dot_thresh <- (MI_dot > 0) & (MI_grat < 0)
	MI_grat_thresh <- (MI_grat > 0) & (MI_dot < 0)

	norm_fits <- list(grat_fits = grat_fits, dot_fits = dot_fits)

	saveRDS(responses_split_list, paste(data_directory, "stimulus_reps.RDS", sep = ""))		
	write.table(epoch_labels, file = paste(data_directory, "epoch_labels.dat", sep = ""), col.names = F, row.names = F)
	saveRDS(norm_fits, paste(data_directory, "response_parametric_fits.RDS", sep=""))

	saveRDS(MI_grat_shuffle, paste(data_directory, "MI_GRAT_shuffle", sep=""))
	saveRDS(MI_dot_shuffle, paste(data_directory, "MI_DOT_shuffle", sep=""))

	write.table(MI_dot_thresh, paste(data_directory, "MI_dot_thresh.dat", sep=""))
	write.table(MI_grat_thresh, paste(data_directory, "MI_grat_thresh.dat", sep=""))

	write.table(MI_dot, paste(data_directory, "MI_grat.dat", sep=""))
	write.table(MI_dot, paste(data_directory, "MI_dot.dat", sep=""))
}

