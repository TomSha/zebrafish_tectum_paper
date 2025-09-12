source("functions/baseline_functions.R")


quant_baseline <- function(prefix){

	data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep = "")
    spikes <- read.table(paste(data_directory, "spikes.dat", sep = ""))
	traces <- np$load(paste(data_directory, "traces_norm.npy", sep = ""))
	norm_fits <- readRDS(paste(data_directory, "response_parametric_fits.RDS", sep=""))

    baseL_index <- calculate_baseline_index(spikes)
    bl_quant <- quant_baseline_metrics(traces, baseL_index)


	write.table(baseL_index, paste(data_directory, "baseL_index.dat", sep = ""), row.names = F, col.names = F)
    write.table(bl_quant$noise, paste(data_directory, "baseline_noise.dat", sep = ""), row.names = F, col.names = F)
    write.table(bl_quant$fluor, paste(data_directory, "baseline_fluor.dat", sep = ""), row.names = F, col.names = F)
}


prefix <- prefix_list[1]
data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep = "")
bl_noise <- read.table(paste(data_directory, "baseline_noise.dat", sep = ""))$V1
MI_dot <- read.table(paste(data_directory, "MI_dot.dat", sep = ""))[,1]
norm_fits <- readRDS(paste(data_directory, "response_parametric_fits.RDS", sep=""))

n_cells <- length(bl_noise)
dot_fits <- norm_fits$dot_fits
signal <- rep(NA, n_cells)

for (i in 1 : n_cells){
    signal[i] <- max(dot_fits[1, i, ])
}

snr <- signal / bl_noise

plot(MI_dot ~ snr)
plot(MI_dot ~ signal)
plot(MI_dot ~ bl_noise)
mod <- lm(MI_dot ~ snr)

response_var <- dot_fits[1,,]
response_mean <- dot_fits[2,,]

plot(log(response_mean) ~ response_var)
a = 4
#plot(log(response_var[,a]) ~ log(response_mean[,a]))
lm(log(response_var[,a]) ~ log(response_mean[,a]))
lm(log(c(response_var)) ~ log(c(response_mean)))




alpha = rep(NA, n_cells)

for (i in 1:n_cells){
    alpha[i] = lm(log(response_var[i,]) ~ log(response_mean[i,]))$coefficients[2]
}


n_stim <- 4

noise <- rexp(n = n_cells)
mu <- rexp(n_cells * n_stim)
sig <- rexp(n_cells * n_stim, rate = 3)

mu <- matrix(mu, nrow = n_cells, ncol = n_stim) + noise

# two separate things
# 1. Are good recording conditions causing MI? - look at noise/baseline
# 2. Are we missing some neurons that have significant mutual information but can#t detect due to noise? SNR?
# because increase in signal can increase SNR and neurons with higher signal are more likely to have higher MI (possibly show with simulation)
# then correlation might be expectee
# 