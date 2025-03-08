library(reticulate)
# Define constants 
np <- import("numpy")
n_planes <- 7	# Number of imaged planes (one fewer than the total number of imaging planes as the 2photon automatically adds a throw away plane)
ac_freq <- 7.28	# The acquisition frequency of the 2photon imaging
delay <- 0 	# The delay between starting imaging and starting the visual stimulation 
n_pixels <- 256	# The n_pixel x n_pixel resolution of the imaging session

# Define subtypes that are included in analysis
inc <- c(2, 3, 4, 5, 7, 22, 8, 23, 10, 11, 12, 14, 1)
groups <- c(1, 2, 3, 4, "5a", "5b", "6a", "6b", "7", 8, 9, 10, 11)
group_names <- c("RC Small", "RC Large", "CR Large", "CR Small", "DS (RC)", "DS (CR)", "SS (Large)", "SS (Small)", "RC SS", "Large DS", "Small DS", "CR SS", "NS")

# Each prefix defines a separate experiment
prefix_list <- read.table("prefix_list")$V1

# The main directory defines the path to where the data from all experiments are stored
# the structure of the data for a given experiment should be: main_directory/prefix/suite2p/combined/
user <- Sys.info()[names(Sys.info()) == "user"]
main_directory <- paste("/media/", user, "/Samsung_T5/Stimulus_Barrage/", sep = "")