source("functions/baseline_functions.R")


quant_baseline <- function(prefix){

	data_directory <- paste(main_directory, prefix, "/suite2p/combined/", sep = "")
    spikes <- read.table(paste(data_directory, "spikes.dat", sep = ""))
	traces <- np$load(paste(data_directory, "traces_norm.npy", sep = ""))

    baseL_index <- calculate_baseline_index(spikes)
    bl_quant <- quant_baseline_metrics(traces, baseL_index)


	write.table(baseL_index, paste(data_directory, "baseL_index.dat", sep = ""), row.names = F, col.names = F)
    write.table(bl_quant$noise, paste(data_directory, "baseline_noise.dat", sep = ""), row.names = F, col.names = F)
    write.table(bl_quant$fluor, paste(data_directory, "baseline_fluor.dat", sep = ""), row.names = F, col.names = F)
}




