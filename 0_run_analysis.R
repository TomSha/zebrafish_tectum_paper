library(reticulate)
library(tiff)
library(parallel)
library(gtools)

source("1_preprocess.R")
source("2_binarise_traces.R")
source("3_calculate_MI.R")

np <- import("numpy")

# Define constants 
n_planes <- 7	# Number of imaged planes (one fewer than the total number of imaging planes as the 2photon automatically adds a throw away plane)
ac_freq <- 7.28	# The acquisition frequency of the 2photon imaging
delay <- 0 	# The delay between starting imaging and starting the visual stimulation 
n_pixels <- 256	# The n_pixel x n_pixel resolution of the imaging session

# Each prefix defines a separate experiment
prefix_list <- read.table("prefix_list")$V1
prefix<-"200122_F1_PVN0" 

# The main directory defines the path to where the data from all experiments are stored
# the structure of the data for a given experiment should be: main_directory/prefix/suite2p/combined/
user <- Sys.info()[names(Sys.info()) == "user"]
main_directory <- paste("/media/", user, "/Samsung_T5/Stimulus_Barrage/", sep = "")

run_analysis <- function(prefix){
	preprocess_dat(prefix)

	binarise_traces(prefix)

	calculate_MI(prefix)
}

