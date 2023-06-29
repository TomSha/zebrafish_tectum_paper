library(partitions)
source("functions/6_neuronal_subtypes_functions.R")



find_subtypes <- function(prefix){
	data_directory <- paste(main_directory,prefix,"/suite2p/combined/",sep="")
	all_responses <- readRDS(paste(data_directory,"stimulus_reps.RDS",sep=""))
	epochs <- read.table(paste(data_directory,"epoch_labels.dat",sep=""))$V1
	xy <- read.table(paste(data_directory,"xy.dat",sep=""))

	ord <- order(epochs)
	all_responses <- all_responses[ord]
	epochs <- epochs[ord]

	# get DOT stimuli
	dot_ep <- grepl("^DOT",epochs)
	dot_responses <- all_responses[dot_ep]
	dot_epochs <- epochs[dot_ep]

	dot_responses <- lapply(dot_responses, function(x) apply(x,c(1,3),mean))
	dot_r <- simplify2array(dot_responses)
	dot_responses <- do.call(cbind,dot_responses)	


	# get GRAT stimuli
	grat_ep <- grepl("GRAT",epochs)
	grat_responses <- all_responses[grat_ep]
	grat_epochs <- epochs[grat_ep]

	grat_responses <- lapply(grat_responses, function(x) apply(x,c(1,3),mean))
	grat_r <- simplify2array(grat_responses)
	grat_responses <- do.call(cbind,grat_responses)	

	n_cell <- nrow(all_responses[[1]])

	dot_clusters <- cluster_cells(dot_responses, dot_epochs)
	grat_clusters <- cluster_cells(grat_responses, grat_epochs)

	dot_model_selec <- rep(NA,n_cell)
	grat_model_selec <- rep(NA,n_cell)

	# for each cell which is best described by spliiting the stimuli into two sets, calculate the selectivity of responses to one set vs the other
	for(i in 1:n_cell){
	       dot_model_selec[i] <- calc_model_selectivity(dot_clusters$model_n[i], dot_r[i, , ], dot_clusters$parti)
	       grat_model_selec[i] <- calc_model_selectivity(grat_clusters$model_n[i], grat_r[i, , ], grat_clusters$parti)
	}

	dot_minus_selec <- dot_model_selec < 0
	dot_minus_selec[is.na(dot_minus_selec)] <- FALSE

	grat_minus_selec <- grat_model_selec < 0
	grat_minus_selec[is.na(grat_minus_selec)] <- FALSE	

	# Any cell which is best described by two subsets, split them into two determined by their response selectivity i.e those that preferentially respond to subset A and those that preferentially respond to subset B
	for(i in dot_clusters$model_2gauss){
		dot_clusters$model_n[dot_clusters$model_n == i & dot_minus_selec] <- i + length(dot_clusters$parti)
	}

	for(i in grat_clusters$model_2gauss){
		grat_clusters$model_n[grat_clusters$model_n == i & grat_minus_selec] <- i + length(grat_clusters$parti)
	}


	write.table(dot_clusters$model_n,paste(data_directory,"dot_model_number.dat",sep=""),col.names=F,row.names=F)
	write.table(grat_clusters$model_n,paste(data_directory,"grat_model_number.dat",sep=""),col.names=F,row.names=F)

	write.table(dot_model_selec,paste(data_directory,"dot_model_selec.dat",sep=""),col.names=F,row.names=F)
	write.table(grat_model_selec,paste(data_directory,"grat_model_selec.dat",sep=""),col.names=F,row.names=F)

	write.table(dot_clusters$LL,paste(data_directory,"dot_LL.dat",sep=""),col.names=F,row.names=F)
	write.table(grat_clusters$LL,paste(data_directory,"grat_LL.dat",sep=""),col.names=F,row.names=F)
}


