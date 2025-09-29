
calculate_baseline_index <- function(spikes, n_bck = 200, n_fwd = 400){

	n_cells = nrow(spikes)
	baseL_index_mat = matrix(nrow = n_cells, ncol = ncol(spikes))

	for (i in 1 : n_cells){
        cat("\r", i)
		spks = spikes[i, ]
		spks_index = which(spks == 1)
		activity_index = lapply(spks_index, function(x) (x - n_bck) : (x + n_fwd))
		activity_index = unique(unlist(activity_index))
		baseL_index = 1 : length(spks)
		baseL_index_mat[i, ] = !baseL_index %in% activity_index
	}

	return(baseL_index_mat)

}

quant_baseline_metrics <- function(traces, baseL_index){

    n_cells = nrow(traces)
    bl_fluor = vector("list", n_cells)
    bl_noise = vector("list", n_cells)
    for (i in 1:n_cells){
        cat("\r", i)
        bl_fluor[[i]] <- mean(traces[i, baseL_index[i, ]])
        bl_noise[[i]] <- sd(traces[i, baseL_index[i, ]])
    }

    bl_fluor <- unlist(bl_fluor)
    bl_noise <- unlist(bl_noise)
    returnList <- list(fluor = bl_fluor, noise = bl_noise)
    
    return(returnList)
}
