
########################################################################
# Rotates 2D coordinates about a point counterclockwise
#
# params:
# alpha = the angle in radians to rotate the points by
# xy = a matrix of the xy coordinates to rotate (neurons x coordinates)
# centre = the central point around which the points are rotated (x,y)
#
# returns:
# xy_r = a matrix of the rotated points (neurons x rotated coordinates)
########################################################################

rotation_mat <- function(alpha, xy, centre = (apply(xy, 2, max) - apply(xy, 2, min))/2 ){

        rotm <- matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), ncol = 2)
        xy_r <- t(rotm %*% (t(xy) - centre) + centre)

	return(xy_r)
}

########################################################################
# Splits up the x axis into n equally spaced bins
# 
# params:
# xy_r = a matrix of the (rotated) xy coordinates of the neurons (neurons x coordinates)
# n_bins = the number of bins (integer)
#
# returns:
# bins = a vector specifying which bin each neuron is in 
########################################################################

split_xaxis <- function(xy_r, n_bins){

	n_cells <- nrow(xy_r)
	xaxis <- xy_r[, 1]
	min_x <- min(xaxis)
	max_x <- max(xaxis)
	bin_w <- seq(from = min_x, to = max_x, length.out = n_bins + 1)
	bins <- rep(NA, n_cells)
	for(i in 1:(n_bins)) bins[xaxis >= bin_w[i] & xaxis < bin_w[i + 1]] <- i
	bins[is.na(bins)] <- n_bins

	return(bins)
}


########################################################################
# takes xy coordinates of neurons and assigns them to one of n equally 
# spaced bins along the anterior-posterior axis of the tectum
# To find the bins we rotate the tectum so that the AP axis is parrallel 
# to the x axis of the image
#
# params:
# xy = a matrix of the xy coordinates (neurons x coordinates)
# angle = the angle in degrees by which the neurons should be rotated
# n_bins = the number of bins to create (integer)
#
# returns:
# bins = a vector specifying which bin each neuron is in 
########################################################################

calculate_bins <- function(xy, angle = 25, n_bins = 10){
	alpha <- angle * pi / 180
	xy_r <- rotation_mat(alpha, xy)
	bins <- split_xaxis(xy_r, n_bins)

	return(bins)
}


########################################################################
# normalises the number of neurons per bin based on the total number 
# of neurons across the bins
#
# params:
# bins = a vector of integers specifying which bin a neuron is in
# thresh = a logical vector which specifies whether a neuron has non-zero mutual information
# NCC_thresh = a logical vector which specifies whether a neuron is visually responsive
# n_bins = the number of bins (we don't find this automatically incase there are no neurons
#          in some bins)
# returns:
# n_cells_norm = a vector of the normalised number of neurons per bin 
########################################################################

norm_cell_num <- function(bins, MI_thresh, NCC_thresh, n_bins = 10){

		n_cells <- table(factor(bins[NCC_thresh], levels = 1:n_bins))
		n_cells_MI <- bins[MI_thresh & NCC_thresh]
		n_cells_MI <- table(factor(n_cells_MI, levels = 1:n_bins))
		n_cells_norm <- n_cells_MI/n_cells
		n_cells_norm[is.na(n_cells_norm)] <- 0

		if(sum(n_cells_norm) == 0){
			return(n_cells_norm)
		}else{
			n_cells_norm <- n_cells_norm / sum(n_cells_norm) * sum(n_cells_MI)
		}
		return(round(n_cells_norm))
}

########################################################################
# Calculates the anterior-posterior bias index for the distribution of 
# neurons across the tectum
#
# params:
# model_output = the output of the hierarchical multinomial model
#
# returns: diff_prob = a list of matrices where each matrix gives the AP bias per experiment (experiment x subtype)
########################################################################

calculate_AP_bias <- function(model_output){

	p <- lapply(model_output, "[[", "p")
	n_stim <- length(model_output)
	stim_names <- names(p)

	# calculate the MAP estimate for each bin along the AP axis for every fish and subtype
	map_p <- lapply(p, function(x) apply(x, c(1, 2, 3), map_estimate))

	# sum the prob. mass for the anterior tectum (bins 1:5) and posterior tectum (bins 6:10) and calc the difference
	diff_prob <- vector("list", n_stim)
	names(diff_prob) <- stim_names

	for(i in 1 : n_stim){
		stim <- map_p[[i]]
		ant_prob <- apply(stim, c(2, 3), function(x) sum(x[1 : 5]))
		post_prob <- apply(stim, c(2, 3), function(x) sum(x[6 : 10]))
		diff_prob[[i]] <- ant_prob - post_prob
	}
	
	return(diff_prob)
}
