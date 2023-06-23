source("constants.R")
source("functions/preprocess_functions.R")

preprocess_dat <- function(prefix){
	
	# Load the data from output of suite2p

	MIP_list <- vector("list", length = n_planes)
	mask_list <- vector("list", length = n_planes)

	data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep="")
	stats_list <- np$load(paste(data_directory, "stat.npy", sep=""), allow_pickle=T)
	traces <- np$load(paste(data_directory, "F.npy", sep=""), allow_pickle=T)
	is_cell_prob <- np$load(paste(data_directory, "/iscell.npy", sep=""), allow_pickle=T)[, 2]
	
	for(i in 1 : n_planes){
		MIP_list[[i]] <- readTIFF(paste(main_directory, prefix, "/", prefix, "_slice", i, "_MIP.tif", sep=""))
		mask_list[[i]] <- readTIFF(paste(main_directory, prefix, "/", prefix, "_slice", i, "_mask.tif", sep=""))
	}

	# extract relevant neurons for analysis 
	# we are finding ROIs that 1) are in the tectum 2) have a non-zero response 3) have prob > 0.5 of being a neuron

	xyz <- get_xyz_coordinates(stats_list, n_pixels)

	tectal_cells_index <- mask_cells(xyz, mask_list)
	responsive_cells_index <- !rowSums(traces) == 0
	is_cells_index <- is_cell_prob > 0.5

	include_index <- is_cells_index * tectal_cells_index * responsive_cells_index

	xyz <- xyz[include_index == 1, ]    	
	traces <- traces[include_index == 1, ]    	

	write.table(xyz, paste(data_directory, "xy.dat", sep=""), row.name=F, col.name=F)
	np$save(paste(data_directory, "traces.npy", sep=""), traces)
}
