library(fitdistrplus)


########################################################################
# Gets the names of the visual stimuli from the timelog file
# 
# params:
# timelog = a data.frame where the 1st column is the time in seconds of when a stimulus started and finished and the 2nd column is the name of the stimulus 
#
# returns:
# epoch_labels = the names of the visual stimuli used
########################################################################

get_epoch_labels <- function(timelog){

	nEP <- nrow(timelog)/2
	epochs <- timelog[, 2][1:nEP*2]
	epoch_labels <- unique(epochs)
	epoch_labels <- epoch_labels[!grepl("CONCENTRIC", epoch_labels)]

	return(epoch_labels)
}

########################################################################
# Splits up the time series into individual epochs of when a stimulus was being shown
# 
# params:
# timelog: a data.frame where the 1st column is the time in seconds of when a stimulus started and finished and the 2nd column is the name of the stimulus 
# traces = a matrix of fluorescence traces (neurons x timepoints)
# ac_freq = imaging acquisition frequency (Hz)
#
# returns:
#
# traces_rep_array_list = a list of 3D arrays where each array is the response to a specific stimulus (neurons x timepoints x repetitions)
# 			  the order of the stimuli is determined by the order of epoch_labels	
########################################################################

split_traces <- function(timelog, traces, ac_freq){

	#split the traces into their epochs

	nEP <- nrow(timelog)/2
	epochs <- timelog[, 2][1:nEP*2]
	epoch_labels <- get_epoch_labels(timelog)		
	timelog[, 1] <- floor((timelog[, 1])*ac_freq)

	start <- timelog[, 1][1:nEP*2-1]
	finish <- timelog[, 1][1:nEP*2]
	
	#make sure each epoch repetition has the same number of frames
	ep_length <- finish-start
	min_ep_length <- sapply(epoch_labels, function(epoch_labels) min(ep_length[grepl(paste("^", epoch_labels, "$", sep = ""), epochs)]))
	for (k in 1:length(epoch_labels)){
			 ep_length[grepl(paste("^", epoch_labels[k], "$", sep = ""), epochs)] <- min_ep_length[k]
	}
	finish <- start + ep_length

	traces_split <- mapply(function(start, finish) traces[1:nrow(traces), start:finish],  start, finish, SIMPLIFY = F)
	traces_split <- lapply(epoch_labels, function(epoch_labels) traces_split[grepl(paste("^", epoch_labels, "$", sep = ""), epochs)])	

	#check each epoch repetition has the same number of frames
	ncols <- sapply(traces_split,  function(x) sapply(x, ncol))

	if(!var(apply(ncols, 2, function(rows) var(rows) == 0)) == 0){
		print("WARNING: the number of frames in your repetitions are not equal, your 3D array will not be correct")
	}
		
	traces_rep_array_list <- lapply(traces_split, simplify2array)

	
	return(traces_rep_array_list)

}

########################################################################
# Fits a univariate normal distribution to the response of each neuron to reptitions of the stimulus
#
# params:
# stimulus_reps_mat = response matrix to repetitions of the same stimulus - the response should be a point estimate of the epoch, eg mean or max response of the cell to that epoch (neurons x responses)
#
# returns:
# est_dist = matrix of mean and std. dev estimates of the response distribution (estimates x neurons)
########################################################################

fit_dist <- function(stimulus_reps_mat){
	est_dist <- apply(stimulus_reps_mat, 1, function(x) fitdist(x, distr = "norm", method = "mme")$estimate)
	return(est_dist)
}



########################################################################
# Calculates the mutual information of a neuron's response to the stimuli as MI(S;R) = H(R) - H(R|S)
# Assuming P(R|S) follows a gaussian distribution
# where S is the stimulus and R is the response of a neuron
#
# params:
# prob =  n_stimuli x 2 matrix which gives the mean (col1) and sd (col2) of a neuron's response to repetitions of a stimulus. 
# gridsize = the size of the grid over which integration is performed
#
# returns:
# MI = the mutual information of the neuron
########################################################################

MI <- function(prob, gridsize = 1000){

	k <- nrow(prob)
	prob <- prob[order(prob[, 1]), ]

	#The range over which the integration is performed (+/- 3sd from smallest to largest mean with n steps)
	x_range <- seq(prob[1, 1] - max(prob[, 2]) * 3, prob[k, 1] + max(prob[, 2]) * 3, length = gridsize)
	channel_matrix <- matrix(NA, k, gridsize)

	#calculate prob.dens. at each x_range value given the prob of response for each k stimulus.
	for(i in 1:k){
		channel_matrix[i, ] <- dnorm(x_range, mean = prob[i, 1], sd = prob[i, 2])
	}

	cond_entropy <- matrix(NA, nrow = k ,ncol = gridsize)
	marg_entropy <- rep(0, gridsize)

	cond_entropy <- (1/k) * channel_matrix * log2(channel_matrix)
	cond_entropy[is.na(cond_entropy)] <- 0
	cond_entropy <- -sum(cond_entropy)
	cond_entropy <- cond_entropy * (x_range[2] - x_range[1])

	channel_matrix_marg <- colSums(channel_matrix) * (1/k)

	marg_entropy <- channel_matrix_marg * log2(channel_matrix_marg)
	marg_entropy[is.na(marg_entropy)] <- 0
	marg_entropy <- -sum(marg_entropy)
	marg_entropy <- marg_entropy * (x_range[2] - x_range[1])

	MI <- marg_entropy - cond_entropy

	return(MI)
}

########################################################################
# Shuffles the responses across stimuli
# 
# params:
# responses_split = a 3D response array (neurons x repetitions x stimuli)
#
# returns:
# responses_shuffle = a list of matrices where each matrix is the shuffled responses to the stimuli (neurons x repetitions)
########################################################################

.shuffle_responses <- function(responses_split){
	# create shuffled responses
	responses_shuffle <- t(apply(responses_split, 1, sample))
	responses_shuffle <- array(responses_shuffle, dim(responses_split))
	
	# create list over the "stimuli" (now shuffled)
	responses_shuffle <- lapply(1:dim(responses_shuffle)[3], function(x) responses_shuffle[, , x])

	return(responses_shuffle)
}

########################################################################
# Calculates the mutual information on shuffled responses.
# This gives us the mutual information assuming the neuron was firing randomly.
#
# params:
# stim_responses_list = a list of matrices where each matrix is the mean response to a specific stimulus (neurons x repetitions)
# n_samples = the number of times the the shuffled MI metric is calculated
# subset_epochs = the stimulus set we want to shuffle (i.e GRAT or DOT)
# epoch_labels = names of the stimuli
#
# returns:
# MI_shuffle_list = a list of vectors where each vector is the MI calculated for every neuron on the shuffled data. The list is length n_samples
########################################################################

calculate_MI_shuffle <- function(stim_responses_list, n_samples = 20){

	stim_responses <- simplify2array(stim_responses_list)
	stim_responses_r <- replicate(stim_responses, n = n_samples, simplify = F)
	MI_shuffle_list <- vector("list", n_samples)

	# shuffle the responses and calculate the MI 
	for(i in 1 : n_samples){
		cat(i,"\r")
		responses_shuffled <- .shuffle_responses(stim_responses_r[[i]])
		norm_fits <- lapply(responses_shuffled, fit_dist)
		norm_fits <- simplify2array(norm_fits)
		MI_shuffle_list[[i]] <- apply(norm_fits, 2, function(x) MI(t(x)))
	}

	return(MI_shuffle_list)
}	
