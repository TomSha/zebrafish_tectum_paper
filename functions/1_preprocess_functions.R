
########################################################################
# define internally called functions

`%inm%` <- function(rows, mat){
	test <- apply(mat, 1, `==`, rows)
  	any(apply(test, 2, all))
}

#checks if a matrix is binary
is_bin_mat <- function(mat){
	is_bin <- identical(as.vector(mat), as.numeric(as.logical(mat)))
	return(is_bin)
}

#convert matrix to binary
convert_2_bin <- function(mat){
	#This assumes the top right pixel should be 0 
	mat[mat == mat[1, n_pixels]] <- 0
	mat[!mat == 0] <- 1
	return(mat)
}
########################################################################

########################################################################
# get xy coordinates function
# This function gets the xyz coordinates of all the neurons
# 
# params:
# stats_list = stat.npy file from suite2p
# n_pixels = the n x n pixel resolution of the image
#
# returns:
# xyz = a n x 3 matrix of xyz coordinates (where n is the number of neurons)
########################################################################

get_xyz_coordinates <- function(stats_list, n_pixels){
	#data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep="")
	#stats_list <- np$load(paste(data_directory, "stat.npy", sep=""), allow_pickle=T)

	xyz <- cbind(t(sapply(stats_list, "[[", "med")), sapply(stats_list, "[[", "iplane"))
	xyz <- matrix(unlist(xyz), ncol = 3)
	xyz[, 3] <- xyz[, 3] + 1

	#for some reason suite2p adds on the number of pixels to some xy coordinates (eg x+256)
	while(sum(xyz > n_pixels) > 0){
		xyz[, 1][xyz[, 1] > n_pixels] <- xyz[, 1][xyz[, 1] > n_pixels] - n_pixels
		xyz[, 2][xyz[, 2] > n_pixels] <- xyz[, 2][xyz[, 2] > n_pixels] - n_pixels
	}

	return(xyz)
}


########################################################################
# mask cells function
# This function finds all neurons from a predetermined mask (e.g those within the tectum)
#
# params:
# xyz = a n x 3 matrix of xyz coordinates
# mask_list = a list where each element is a n_pixel x n_pixel binary mask of the tectum for each imaging plane
#
#returns:
#tectal_cell_index = a logical vector specifying whether a neuron is in the tectum or not
########################################################################

mask_cells <- function(xyz, mask_list){   
	
	n_planes <- length(mask_list)
	tectal_cells_index_list <- vector("list", n_planes + 1)

	for (i in 1 : n_planes){
		mask <- mask_list[[i]]
		
		if(!is_bin_mat(mask)){
			print("WARNING: mask is not binary,  converting to binary")
			mask <- convert_2_bin(mask)
		}
		x <- floor(xyz[, 1])
		x <- x[xyz[, 3] == i]
		x[x == 0] <- 1
		y <- floor(xyz[, 2])
		y[y == 0] <- 1
		y <- y[xyz[, 3] == i]

		cells_vec <- rep(0, n_pixels * n_pixels)
		cells_vec[x + n_pixels * (y - 1)] <- 1
		cells_mat <- matrix(cells_vec, ncol = n_pixels)
		tectal_cells <- mask * cells_mat 
		tectal_cells <- which(tectal_cells == 1, arr.ind = T)
		tectal_cells_index_list[[i]] <- apply(cbind(x, y), 1, `%inm%`, tectal_cells)		

}
	#we need to set the cells in the throw away imaging plane to FALSE in the tectal cell index
	tectal_cells_index_list[[i + 1]] <- rep(F, sum(xyz[, 3] == n_planes + 1))
	tectal_cells_index <- unlist(tectal_cells_index_list)
	
	return(tectal_cells_index)
}





