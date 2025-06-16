source("constants.R")
source("functions/3_MI_functions.R")
source("functions/10_decoding_functions.R")
source("functions/4_ant_post_bias_functions.R")
library(gam)
library(tiff)

calc_LOO_logliklihood <- function(prefix_list){

	for (i in 1 : length(prefix_list)){

		prefix <- prefix_list[i]
		print(prefix)

		data_directory <- paste(main_directory,prefix,"/suite2p/combined/",sep="")
		timelog <- read.table(paste(main_directory,prefix,"/time_epoches.log",sep=""))
		responses_split_list <- readRDS(paste(data_directory, "stimulus_reps.RDS", sep = ""))
		
		responses_mean_list <- lapply(responses_split_list, function(x) apply(x, c(1, 3), mean))
		epoch_labels <- get_epoch_labels(timelog)
		responses_list <- responses_mean_list[grep("GRAT|^DOT", epoch_labels)]
		n_stim <- length(responses_list)

		fits <- lapply(responses_list, LOO_fit)
		fits <- simplify2array(fits)
		responses_array <- simplify2array(responses_list)

		log_lik <- calc_likelihood(responses_array, fits)
		saveRDS(log_lik, paste(data_directory, "LOO_loglik.RDS", sep = ""))
		
	}
}


for (i in 1:length(prefix_list)){
	pvn_cells_index_list <- vector("list", 9)
	prefix <- prefix_list[i]
	data_directory <- paste(main_directory,prefix,"/suite2p/combined/",sep="")
	xy <- read.table(paste(data_directory, "xy.dat", sep = ""))
	xy_r <- read.table(paste(main_directory,prefix,"/xy_warped_2.dat",sep=""))
	for (j in 1:7){
		mask <- readTIFF(paste(main_directory, "/slice_",j,"_PVN_mask.tif", sep=""))
		mask <- rotate(mask)
		pvn_cells_index_list[[j]] <- mask_cells(xy_r[xy[,3] == j, 1:2], mask)
	}
	pvn_cells_index <- unlist(pvn_cells_index_list)
	write.table(pvn_cells_index, paste(data_directory, "pvn_cells_index.dat", sep = ""), col.names = F, row.names = F)

}


for (prefix in prefix_list){
	print(prefix)
	data_directory <- paste(main_directory,prefix,"/suite2p/combined/",sep="")
	log_lik <- readRDS(paste(data_directory, "LOO_loglik.RDS", sep = ""))
	timelog <- read.table(paste(main_directory,prefix,"/time_epoches.log",sep=""))
	xy <- read.table(paste(data_directory, "xy.dat", sep = ""))
	pvn_cells_index <- read.table(paste(data_directory, "pvn_cells_index.dat", sep = ""))$V1
	NCC_thresh <- read.table(paste(data_directory, "NCC_thresh.dat", sep = ""))[,1]

	xy <- xy[pvn_cells_index & NCC_thresh,]
	log_lik <- log_lik[,,pvn_cells_index & NCC_thresh]
#	xy <- xy[pvn_cells_index, ]
#	log_lik <- log_lik[,,pvn_cells_index]


	epoch_labels <- get_epoch_labels(timelog)
	epoch_labels <- epoch_labels[grep("GRAT|^DOT", epoch_labels)]

	dot_labels <- grep("^DOT", epoch_labels)
	grat_labels <- grep("GRAT", epoch_labels)

	get_indices <- function(x) (((x - 1) * 10) + 1) : (((x - 1) * 10) + 10)
	dot_indices <- c(sapply(dot_labels, get_indices))
	grat_indices <- c(sapply(grat_labels, get_indices))

	dot_log_lik <- log_lik[dot_indices, dot_labels,]
	grat_log_lik <- log_lik[grat_indices, grat_labels, ]
	k <- 100 
#	ks <- c(1 : 20, seq(30, k, 10))
	ks <- 1 : k
	dot_decoded <- decode_stimulus(dot_log_lik, xy, k, ks)
	grat_decoded <- decode_stimulus(grat_log_lik, xy, k, ks)

	saveRDS(dot_decoded, paste(data_directory, "ML_dot_decoded.RDS", sep = ""))
	saveRDS(grat_decoded, paste(data_directory, "ML_grat_decoded.RDS", sep = ""))
	
}




NCC_thresh_list <- vector("list", length(prefix_list))
min_max_x_list <- vector("list", length(prefix_list))
xy_r_list <- vector("list", length(prefix_list))
xy_list <- vector("list", length(prefix_list))
pvn_cells_index_list <- vector("list", length(prefix_list))
xy_r_projection_list <- vector("list", length(prefix_list))
dot_decoded_correct_list <- vector("list", length(prefix_list))
dot_stim_decoded_list <- vector("list", length(prefix_list))
grat_decoded_correct_list <- vector("list", length(prefix_list))
grat_stim_decoded_list <- vector("list", length(prefix_list))
dot_average_list <- vector("list", length(prefix_list))
grat_average_list <- vector("list", length(prefix_list))
dot_labels_list <- vector("list", length(prefix_list))
grat_labels_list <- vector("list", length(prefix_list))

angle <- 25
alpha <- angle * pi / 180

n_reps <- 10
actual <- rep(1: 4, each = n_reps)

for (i in 1:length(prefix_list)){
	prefix <- prefix_list[i]
	data_directory <- paste(main_directory,prefix,"/suite2p/combined/",sep="")
	NCC_thresh_list[[i]] <- read.table(paste(data_directory, "NCC_thresh.dat", sep = ""))[,1]
	pvn_cells_index_list[[i]] <- read.table(paste(data_directory, "pvn_cells_index.dat", sep = ""))$V1	
	xy_r_list[[i]] <- read.table(paste(main_directory,prefix,"/xy_warped_2.dat",sep=""))
	xy_list[[i]] <- read.table(paste(data_directory, "xy.dat", sep = ""))

	xy_r_projection_list[[i]] <- rotation_mat(alpha,xy_r_list[[i]][, 1 : 2], alpha)
	min_max_x_list[[i]] <- c(min(xy_r_projection_list[[i]][, 1]), max(xy_r_projection_list[[i]][, 1]))

	xy_r_projection_list[[i]] <- xy_r_projection_list[[i]][pvn_cells_index_list[[i]] & NCC_thresh_list[[i]],]
	xy_r_list[[i]] <- xy_r_list[[i]][pvn_cells_index_list[[i]] & NCC_thresh_list[[i]],]
	xy_list[[i]] <- xy_list[[i]][pvn_cells_index_list[[i]] & NCC_thresh_list[[i]],]

	dot_decoded <- readRDS(paste(data_directory, "ML_dot_decoded.RDS", sep = ""))
	grat_decoded <- readRDS(paste(data_directory, "ML_grat_decoded.RDS", sep = ""))

	dot_decoded_correct_list[[i]] <- apply(dot_decoded, c(1, 3), function(x) sum(x == actual) / length(x))
	grat_decoded_correct_list[[i]] <- apply(grat_decoded, c(1, 3), function(x) sum(x == actual) / length(x))
	dot_stim_decoded_list[[i]] <- calc_stim_decoding(actual, dot_decoded)
	grat_stim_decoded_list[[i]] <- calc_stim_decoding(actual, grat_decoded)
	timelog <- read.table(paste(main_directory,prefix,"/time_epoches.log",sep=""))

	epoch_labels <- get_epoch_labels(timelog)

	dot_labels_list[[i]] <- epoch_labels[grep("^DOT", epoch_labels)]
	grat_labels_list[[i]] <- epoch_labels[grep("GRAT", epoch_labels)]

	dot_ord <- order(dot_labels_list[[i]])
	grat_ord <- order(grat_labels_list[[i]])

	dot_stim_decoded_list[[i]] <- dot_stim_decoded_list[[i]][dot_ord]
	grat_stim_decoded_list[[i]] <- grat_stim_decoded_list[[i]][grat_ord]
	dot_labels_list[[i]] <- dot_labels_list[[i]][dot_ord]
	grat_labels_list[[i]] <- grat_labels_list[[i]][grat_ord]

}

		


plot_decoded_correct_perc <- function(S_lab = 1.5, S_axis = 1.5, S_leg = 2, L1 = 3, L2 = 3, ks = 1:100){

	mu_grats <- sapply(grat_decoded_correct_list, function(x) apply(x, 1, mean))
	mu_dots <- sapply(dot_decoded_correct_list, function(x) apply(x, 1, mean))
	matplot(ks, mu_dots * 100, col = add_alpha(0.5, martincolourscale[1]), 
	pch = 19, ylab = "", ylim = c(0, 100), xlab = "", font.axis = 2, cex.axis = S_axis)
	matlines(ks, mu_dots * 100, col = martincolourscale[1], lty = 1)

	matplot(ks, mu_grats * 100, col = add_alpha(0.5, martincolourscale[2]), pch = 19, add = TRUE)
	matlines(ks, mu_grats * 100, col = martincolourscale[2], lty = 1)
	legend("bottomright", fill = martincolourscale[1:2], legend = c("local motion", "whole-field motoin"), bty = "n", cex = S_leg, text.font = 2)

	mtext(side = 1, text = "neurons/assembly", cex = S_lab, font = 2, line = L1)
	mtext(side = 2, text = "decoding performance (%)", cex = S_lab, font = 2, line = L2)

}

plot_ensemble <- function(xy_r, NN, neuron = 7){
	plot(xy_r[, 1:2], pch = 19, col = rgb(0, 0, 0, 0.1), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
	points(xy_r[NN[, neuron], 1:2], pch = 19)

}


fit_gam <- function(x, y, min_x = NA, max_x = NA, cols = martincolourscale[1], plot_gam = TRUE){

		if(is.na(min_x)) min_x <- min(x)
		if(is.na(max_x)) max_x <- max(x)

		fit <- gam(y ~ s(x))

		new_dat <- data.frame(x = min_x:max_x)
		preds <- predict(fit, new_dat)
#		ord <- order(x)
		if(plot_gam){
			points(new_dat[, 1], (preds * 100), type = "l", col = cols, lwd = 3, ylab = "", xlab = "")
#			abline(v = new_dat[which.max(preds),1], lty = 2, lwd = 2, col = add_alpha(1, cols))
		}else{
			return(preds)
		}

}


plot_average_regression <- function(xy_r_projection_list, dot_decoded_correct_list, grat_decoded_correct_list, k,
S_axis = 1.5, S_main = 1.5, S_lab = 1.5, S_leg = 2, L1 = 2, L2 = 3, min_x = -70, max_x = 150){

	dot_preds_mat <- matrix(nrow = length(prefix_list), ncol = length(min_x : max_x))
	grat_preds_mat <- matrix(nrow = length(prefix_list), ncol = length(min_x : max_x))

	for (i in 1:length(prefix_list)){
		dot_preds_mat[i, ] <- fit_gam(xy_r_projection_list[[i]][,1], dot_decoded_correct_list[[i]][k,], plot_gam = FALSE, min_x = min_x, max_x = max_x)
		grat_preds_mat[i, ] <- fit_gam(xy_r_projection_list[[i]][,1], grat_decoded_correct_list[[i]][k,], plot_gam = FALSE, min_x = min_x, max_x = max_x)
	}

	mu <- colMeans(dot_preds_mat * 100)
	sig <- apply(dot_preds_mat * 100, 2, sd)
	upper <- mu + sig
	lower <- mu - sig

	plot(min_x : max_x, mu, col = martincolourscale[1], type = "l", lwd = 2,
	ylim = c(0, 100), xlim = c(min_x, max_x), xaxt = "n", cex.axis = S_axis, cex.main = S_main, font.axis = 2, ann = FALSE)
	polygon(c(min_x:max_x, max_x:min_x), c(upper, rev(lower)), col = add_alpha(0.3, martincolourscale[1]), border = FALSE)
#	abline(v = (min_x : max_x)[which.max(mu)], lty = 2, lwd = 2, col = martincolourscale[1])


	mu <- colMeans(grat_preds_mat * 100)
	sig <- apply(grat_preds_mat * 100, 2, sd)
	upper <- mu + sig
	lower <- mu - sig

	points(min_x : max_x, mu, col = martincolourscale[2], type = "l", lwd = 2)
	polygon(c(min_x:max_x, max_x:min_x), c(upper, rev(lower)), col = add_alpha(0.3, martincolourscale[2]), border = FALSE)
#	abline(v = (min_x : max_x)[which.max(mu)], lty = 2, lwd = 2, col = martincolourscale[2])

#	abline(v = (max_x - abs(min_x))/2, lty = 2, lwd = 2, col = rgb(0, 0, 0, 0.3))

	mtext(side = 1, text = "anterior:posterior", cex = S_lab, font = 2, line = L1)
	mtext(side = 2, text = "decoding performance (%)", cex = S_lab, font = 2, line = L2)

	legend("bottomright", fill = martincolourscale[1:2], legend = c("local motion", "whole-field motoin"), bty = "n", cex = S_leg, text.font = 2)

}




plot_location_v_performance <- function(xy_r_projection, dot_decoded_correct, grat_decoded_correct, k, reg_line = TRUE,
S_lab = 1.5, L2 = 3, L1 = 2, S_main = 3, S_axis = 1.5, main_t = "", min_x = -70, max_x = 150){

	xy_r_projection <- xy_r_projection[,1]
	plot(xy_r_projection, dot_decoded_correct[k,] * 100, cex.axis = S_axis, cex.lab = S_lab, cex.main = S_main,
	pch = 19, ylim = c(0, 100), col = add_alpha(0.2, martincolourscale[1]), xaxt = "n", ylab = "", 
	xlab = "", font.axis = 2, main = main_t, xlim = c(min_x, max_x))

	points(xy_r_projection, grat_decoded_correct[k,] * 100, 
	pch = 19, ylim = c(0, 100), col = add_alpha(0.2, martincolourscale[2]), xaxt = "n")

	if(reg_line){
		fit_gam(xy_r_projection, dot_decoded_correct[k,], cols = martincolourscale[1])
		fit_gam(xy_r_projection, grat_decoded_correct[k,], cols = martincolourscale[2])
	}

#	abline(v = max_x/2, lty = 2, lwd = 2, col = rgb(0, 0, 0, 0.3))

	mtext(side = 1, text = "anterior:posterior", cex = S_lab, font = 2, line = L1)
	mtext(side = 2, text = "decoding performance (%)", cex = S_lab, font = 2, line = L2)

}

calc_bias <- function(decoded_correct, xy_r_projection, min_max_x){
	min_x <- min_max_x[1]
	max_x <- min_max_x[2]
	half <- (max_x - abs(min_x)) / 2

	xy_r_projection <- xy_r_projection[,1]

	decoded_correct <- decoded_correct * 100

	ant_mean <- rep(NA, nrow(decoded_correct))
	post_mean <- rep(NA, nrow(decoded_correct))

	for (i in 1 : nrow(decoded_correct)){
		post_mean[i] <-	mean(decoded_correct[i, xy_r_projection > half])
		ant_mean[i] <-	mean(decoded_correct[i, xy_r_projection <= half])
	}
	 return(ant_mean - post_mean)
	
}

plot_bias <- function(xy_r_projection_list, dot_decoded_correct_list, grat_decoded_correct_list, S_lab = 1.5, S_axis = 1.5, 
L2 = 3, L1 = 3, S_leg = 2, ylims = c(-25, 25), main_t = "", S_main = 3, ks = 1:100){
	dot_ant_post_bias <- mapply(calc_bias, dot_decoded_correct_list, xy_r_projection_list, min_max_x_list)
	grat_ant_post_bias <- mapply(calc_bias, grat_decoded_correct_list, xy_r_projection_list, min_max_x_list)

	matplot(ks, dot_ant_post_bias, pch = 19, col = add_alpha(0.3, martincolourscale[1]), ylim = ylims,
	xlab = "", ylab = "", font.axis = 2, cex.axis = S_axis, main = main_t, cex.main = S_main)
	matlines(ks, dot_ant_post_bias, lty = 1, pch = 19, col = add_alpha(0.3, martincolourscale[1]))
	points(ks, rowMeans(dot_ant_post_bias), pch = 19, col = martincolourscale[1], type = "l", lwd = 3)

	matplot(ks, grat_ant_post_bias, pch = 19, col = add_alpha(0.3, martincolourscale[2]), add = T)
	matlines(ks, grat_ant_post_bias, lty =1,  pch = 19, col = add_alpha(0.3, martincolourscale[2]))
	points(ks, rowMeans(grat_ant_post_bias), pch = 19, col = martincolourscale[2], type = "l", lwd = 3)
	abline(h = 0, lty = 2, lwd = 2, col = rgb(0, 0, 0, 0.3))

	mtext(side = 1, text = "neurons/assembly", cex = S_lab, font = 2, line = L1)
	mtext(side = 2, text = "anterior:posterior bias (%)", cex = S_lab, font = 2, line = L2)

	legend("bottomright", fill = martincolourscale[1:2], legend = c("local motion", "whole-field motoin"), bty = "n", cex = S_leg, text.font = 2)
}



plot_fig <- function(fish = 7, neuron = 1, S_letter = 3, S_lab = 2, S_axis = 2){
	
	par(mfrow = c(3, 3))
	NN <- find_KNN(xy_list[[fish]], 200)
	par(mar = c(3, 6, 3, 6))
	plot_ensemble(xy_r_list[[fish]], NN, neuron)
	title(main="A)",adj=0,cex.main=S_letter)
	par(mar = c(6, 6, 5, 2))
	plot_decoded_correct_perc(S_lab = S_lab, S_axis = S_axis)
	title(main="B)",adj=0,cex.main=S_letter)
#	plot(xy_r_list[[fish]][, 1:2], col = colPalette(xy_r_projection_list[[fish]]), pch = 19, ann = FALSE, xaxt = "n", yaxt = "n")
	plot_bias(xy_r_projection_list, dot_decoded_correct_list, grat_decoded_correct_list, S_lab = S_lab, S_axis = S_axis)
	title(main="C)",adj=0,cex.main=S_letter)

	plot_location_v_performance(xy_r_projection_list[[fish]], dot_decoded_correct_list[[fish]], grat_decoded_correct_list[[fish]], k = 1, main_t = "single neuron decoding", S_lab = S_lab, S_axis = S_axis)
	title(main="D)",adj=0,cex.main=S_letter)
	plot_location_v_performance(xy_r_projection_list[[fish]], dot_decoded_correct_list[[fish]], grat_decoded_correct_list[[fish]], k = 10, main_t = "10 neurons/assembly", S_lab = S_lab, S_axis = S_axis)
	title(main="E)",adj=0,cex.main=S_letter)
	plot_location_v_performance(xy_r_projection_list[[fish]], dot_decoded_correct_list[[fish]], grat_decoded_correct_list[[fish]], k = 20, main_t = "20 neurons/assembly", S_lab = S_lab, S_axis = S_axis)
	title(main="F)",adj=0,cex.main=S_letter)

	plot_average_regression(xy_r_projection_list, dot_decoded_correct_list, grat_decoded_correct_list, 1, S_lab = S_lab, S_axis = S_axis)
	title(main="G)",adj=0,cex.main=S_letter)
	plot_average_regression(xy_r_projection_list, dot_decoded_correct_list, grat_decoded_correct_list, 2, S_lab = S_lab, S_axis = S_axis)
	title(main="H)",adj=0,cex.main=S_letter)
	plot_average_regression(xy_r_projection_list, dot_decoded_correct_list, grat_decoded_correct_list, 40, S_lab = S_lab, S_axis = S_axis)
	title(main="I)",adj=0,cex.main=S_letter)

}



plot_fig_stim <- function(S_letter = 3, S_main = 3, S_lab = 2, S_axis = 2){

	par(mfrow=c(1,4), mar = c(5, 5, 5, 2))
	#check the order!!!!!
	stim_names <- c("CR 5", "CR 26", "RC 5", "RC 26")
	fig_letters <- c("J)", "K)", "L)", "M)")

	for(i in 1:4){
		plot_bias(xy_r_projection_list, lapply(dot_stim_decoded_list, "[[", i), lapply(grat_stim_decoded_list, "[[", i), ylims = c(-35, 35), main_t = paste("stimulus:", stim_names[i]), S_main = S_main, S_lab = S_lab, S_axis = S_axis)
		title(main=fig_letters[i],adj=0,cex.main=S_letter)
	}
}



user<-Sys.info()[names(Sys.info())=="user"]
save_direc<-paste("/media/",user,"/Samsung_T5/Work/tectal_subdomain_paper/",sep="")
png(paste(save_direc, "fig_10.png", sep = ""), width = 1500, height = 1000)
plot_fig(2,20)
dev.off()

png(paste(save_direc, "fig_10b.png", sep = ""), width = 1500, height = 1000/3)
plot_fig_stim()
dev.off()




######## validation
# response of neuron 1 for repetition i to stimulus n
#neuron = 560
#repe = 5
#stimulus = 7
#
#r <- responses_array[neuron, repe, stimulus]
#
## params of neuron 1 to repetition i and stimulus n
#mu <- fits[1, neuron, repe,]
#sig <- fits[2, neuron, repe,]
#
#dnorm(r, mu, sig, log = TRUE)
#log_lik[stimulus * 10 - 10 + repe, , neuron]
#log_lik[20,]
#


plot_halfway <- function(xy_r_projection, min_max_x){
	half <- (max_x - abs(min_x)) / 2

	plot(xy_r_projection, pch = 19, col = "white")
	points(xy_r_projection[xy_r_projection[,1] > half, ], pch = 19, col = "blue")
	points(xy_r_projection[xy_r_projection[,1] <= half, ], pch = 19, col = "red")	
}

par(mfrow = c(3,3))
mapply(function(xy_r, min_max) plot_halfway(xy_r, min_max), xy_r_projection_list, min_max_x_list)