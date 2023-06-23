source("functions/2_binarise_traces_functions.R")

# The function takes the time series for all tectal neurons and binarises them.
# We then find all neurons which are visually responsive.
# Visually responsive neurons are defined as those neurons whose firing correlates with the presentation of the visual stimuli.

binarise_traces <- function(prefix){

	data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep = "")
	traces <- np$load(paste(data_directory, "traces.npy", sep = ""), allow_pickle = T)
	timelog <- read.table(paste(main_directory, prefix, "/time_epoches.log", sep = ""))

	ncells <- nrow(traces)
	nz <- ncol(traces)

	bin_traces <- matrix(nrow = ncells, ncol = nz)
	traces <- traces / rowMeans(traces)

	# Here we binarise the traces using bcl function
	bb <- lapply(1 : ncells,function(x) bcl(traces[x,], 0.01, 0.8, 9e-3, 1e-3, .4, 0, 20, "EXP"))
	spikes <- t(sapply(bb, "[[", 3))
	bin_baseline <- t(sapply(bb, "[[", 2))

	# Here we subtract baseline drift from the traces
	traces <- traces - bin_baseline

	# Here we use the binarised traces to calculate the correlation coefficient between the firing of each neuron and the presentation of the stimuli
	corr <- calculate_correlation(timelog, spikes, delay, ac_freq)

	np$save(paste(main_directory,prefix,"/suite2p/combined/traces_norm.npy",sep=""),traces)	
	write.table(corr$corvec, paste(data_directory, "corvec.dat", sep = ""), row.names = F, col.names = F)
	write.table(corr$threshMean, paste(data_directory, "corvec_thresh.dat", sep = ""), row.names = F, col.names = F)
	write.table(corr$threshSD, paste(data_directory, "corvec_threshSD.dat", sep = ""), row.names = F, col.names = F)
	write.table(corr$thresh, paste(data_directory, "NCC_thresh.dat", sep = ""), row.names = F, col.names = F)
}
