library(reticulate)
library(tiff)
library(parallel)
library(gtools)

source("1_preprocess.R")

# Define constants 
n_planes <- 7
ac_freq <- 7.28
delay <- 0
n_pixels <- 256

# read in the prefix which defines each experiment
prefix_list <- read.table("prefix_list")$V1
prefix<-"200122_F1_PVN0" 

np <- import("numpy")
user <- Sys.info()[names(Sys.info()) == "user"]
main_directory <- paste("/media/", user, "/Samsung_T5/Stimulus_Barrage/", sep = "")

preprocess_dat(prefix)



