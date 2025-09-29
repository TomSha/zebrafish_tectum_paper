

MI_dot_list<-vector("list",length(prefix_list))
MI_grat_list<-vector("list",length(prefix_list))

MI_dot_thresh_list<-vector("list",length(prefix_list))
MI_grat_thresh_list<-vector("list",length(prefix_list))

NCC_thresh_list<-vector("list",length(prefix_list))

baseline_noise_list <- vector("list", length(prefix_list))
baseline_fluor_list <- vector("list", length(prefix_list))
snr_list <- vector("list", length(prefix_list))

for(i in 1:length(prefix_list)){
	print(i)
	prefix<-prefix_list[i]
	data_directory<-paste(main_directory,prefix,"/suite2p/combined/",sep="")

	MI_dot_list[[i]]<-read.table(paste(data_directory,"MI_dot.dat",sep=""))[,1]
	MI_grat_list[[i]]<-read.table(paste(data_directory,"MI_grat.dat",sep=""))[,1]

	MI_grat_thresh_list[[i]]<-read.table(paste(data_directory,"MI_grat_thresh.dat",sep=""))[,1]
	MI_dot_thresh_list[[i]]<-read.table(paste(data_directory,"MI_dot_thresh.dat",sep=""))[,1]

	NCC_thresh_list[[i]]<-read.table(paste(data_directory,"NCC_thresh.dat",sep=""))[,1]

    baseline_noise_list[[i]]<-read.table(paste(data_directory,"baseline_noise.dat",sep=""))[,1]
    baseline_fluor_list[[i]]<-read.table(paste(data_directory,"baseline_fluor.dat",sep=""))[,1]
    snr_list[[i]]<-baseline_fluor_list[[i]] / baseline_noise_list[[i]]

}





plot_noise_vs_dot_MI <- function(i, thresh = TRUE){

    if(thresh){
        MI <- MI_dot_list[[i]][MI_dot_thresh_list[[i]] & NCC_thresh_list[[i]]]
        noise <- baseline_noise_list[[i]][MI_dot_thresh_list[[i]] & NCC_thresh_list[[i]]]
    }else{
        MI <- MI_dot_list[[i]]
        fluor <- baseline_noise_list[[i]]
    }

   plot(noise, MI, pch = 19, col = rgb(0, 0, 0, 0.3)) 

}

plot_fluor_vs_dot_MI <- function(i, thresh = TRUE){

    if(thresh){
        MI <- MI_dot_list[[i]][MI_dot_thresh_list[[i]] & NCC_thresh_list[[i]]]
        fluor <- baseline_fluor_list[[i]][MI_dot_thresh_list[[i]] & NCC_thresh_list[[i]]]
    }else{
        MI <- MI_dot_list[[i]]
        fluor <- baseline_fluor_list[[i]]
    }
   plot(fluor, MI, pch = 19, col = rgb(0, 0, 0, 0.3)) 
   }

plot_snr_vs_dot_MI <- function(i, thresh = TRUE){

    if(thresh){
        MI <- MI_dot_list[[i]][MI_dot_thresh_list[[i]] & NCC_thresh_list[[i]]]
        snr <- snr_list[[i]][MI_dot_thresh_list[[i]] & NCC_thresh_list[[i]]]
    }else{
        MI <- MI_dot_list[[i]]
        snr <- snr_list[[i]]
    }
   plot(MI ~ snr, pch = 19, col = rgb(0, 0, 0, 0.3)) 
   }

plot_noise_vs_grat_MI <- function(i, thresh = TRUE){

    if(thresh){
        MI <- MI_grat_list[[i]][MI_grat_thresh_list[[i]] & NCC_thresh_list[[i]]]
        noise <- baseline_noise_list[[i]][MI_grat_thresh_list[[i]] & NCC_thresh_list[[i]]]
    }else{
        MI <- MI_grat_list[[i]]
        fluor <- baseline_noise_list[[i]]
    }

   plot(noise, MI, pch = 19, col = rgb(0, 0, 0, 0.3)) 

}

plot_fluor_vs_grat_MI <- function(i, thresh = TRUE){

    if(thresh){
        MI <- MI_grat_list[[i]][MI_grat_thresh_list[[i]] & NCC_thresh_list[[i]]]
        fluor <- baseline_fluor_list[[i]][MI_grat_thresh_list[[i]] & NCC_thresh_list[[i]]]
    }else{
        MI <- MI_grat_list[[i]]
        fluor <- baseline_fluor_list[[i]]
    }
   plot(fluor, MI, pch = 19, col = rgb(0, 0, 0, 0.3)) 
   }


plot_basline_func <- function(i){
    par(mfrow = c(2,1))
    plot(traces[i,], type = "l")
    lines(which(spikes[i,] == 1), rep(50, length(which(spikes[i,] ==1))), type = "p", col = "red", pch = 19)
    if (sum(baseL_index[i,]) != 0){
        lines(which(baseL_index[i,]), traces[i,baseL_index[i,]], col = "green4")
        plot(traces[i, baseL_index[i,]], type = "l") 
    }


}