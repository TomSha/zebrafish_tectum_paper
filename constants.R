library(reticulate)
library(tiff)
library(parallel)
library(gtools)

source("~/2photon/plot_functions/martincolourscale.R")

n_planes <- 7
ac_freq <- 7.28
delay <- 0
n_pixels <- 256

rotate <- function(x) t(apply(x, 2, rev))


np<-import("numpy")
user<-Sys.info()[names(Sys.info())=="user"]
main_directory<-paste("/media/",user,"/Samsung_T5/Stimulus_Barrage/",sep="")

pulled_directory<-paste(main_directory,"pulled/",sep="")

titles<-c("CR 5","RC 5","5","CR 26","RC 26","26","RC","CR")
elements<-c(5,2,23,4,3,8,7,22)

inc <- c(2, 3, 4, 5, 7, 22, 8, 23, 10, 11, 12, 14, 1)

groups <- c(1, 2, 3, 4, "5a", "5b", "6a", "6b", "7", 8, 9, 10, 11)

group_names <- c("RC Small", "RC Large", "CR Large", "CR Small", "DS (RC)", "DS (CR)", "SS (Large)", "SS (Small)", "RC SS", "Large DS", "Small DS", "CR SS", "NS")

empty_plot<-function(){
	plot(1, type="n", xlab="", ylab="",bty="n",xaxt="n",yaxt="n")
	axis(1,col=rgb(0, 0, 0, 0), col.axis=rgb(0, 0, 0, 0))
}
