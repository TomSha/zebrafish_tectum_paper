library(reticulate)
library(tiff)
library(parallel)
library(gtools)

source("1_preprocess.R")
source("2_binarise_traces.R")
source("3_calculate_MI.R")
source("4_ant_post_bias.R")
source("5_gaussian_process.R")
source("6_neuronal_subtypes.R")
source("7_ant_post_bias_subtypes.R")
source("8_linear_regression.R")
source("9_quantify_baseline.R")
source("constants.R")



run_analysis <- function(prefix_list){

	for(i in 1:prefix_list){
		prefix <- prefix_list[i]
		preprocess_dat(prefix)	# 1

		binarise_traces(prefix) # 2

		calculate_MI(prefix) 	# 3
	}

	AP_bias(prefix_list)		# 4

	for(i in 1 : prefix_list){
		prefix <- prefix_list[i]

		save_GP_fit(prefix)	# 5

		find_subtypes(prefix)	# 6
	}
		
	AP_bias_subtypes(prefix_list)	# 7
	save_hier_lin_reg(prefix_list)	# 8

	for(i in 1 : prefix_list){

		quant_baseline(prefix)
	}

}

