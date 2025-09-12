# estimate P(R|S) for each neuron and stimulus leaving out 1 repetition
LOO_fit <- function(stimulus_reps_mat){

	n_reps <- ncol(stimulus_reps_mat)
	est_dist_list <- vector("list", n_reps)
    print("LOO fit")

  for (i in 1 : n_reps){
    cat(i,"\r")
    
    est_dist_list[[i]] <- fit_dist(stimulus_reps_mat[, -i])
    }
	
	return(simplify2array(est_dist_list))

}

fit_multi_norm <- function(dat){

		mu_hat <- colMeans(dat)
		cor_mat <- cor(dat)
		var_vec <- apply(dat, 2, sd)
		var_mat <- matrix(rep(var_vec, each = length(mu_hat)), nrow = length(mu_hat))
		Sigma_hat <- cor_mat * var_mat * t(var_mat)

		returnList <- list("mu_hat" = mu_hat, "Sigma_hat" = Sigma_hat)
		
		return(returnList)
}

LOO_multi_fit <- function(dat){

		params_list <- vector("list", length(n_reps))
		for (i in 1 : n_reps){
			params_list[[i]] <- fit_multi_norm(t(dat[,-i]))
		}
		return(params_list)
}


calc_likelihood <- function(responses_array, fits){

    n_cells <- dim(responses_array)[1]
    n_stim <- dim(responses_array)[3]

    log_lik_list <- vector("list", n_cells)
    print("calc likelihood")

    for (i in 1 : n_cells){
        cat(i,"\r")
        
        # response of neuron i for each repetition x stimulus
        r <- responses_array[i, , ]
        r <- rep(c(r), each = n_stim)
        mu <- fits[1, i, , ]
        sig <- fits[2, i, , ]

        mu <- rep(c(t(mu)), n_stim)
        sig <- rep(c(t(sig)), n_stim)

        log_lik <- dnorm(r, mu, sig, log = TRUE)
        log_lik_list[[i]] <- matrix(log_lik, ncol = n_stim, byrow = TRUE)
    }

    log_lik <- simplify2array(log_lik_list)
    return(log_lik)
}


# find the k nearest neighbours of all cell
# xy = xy (and z) coordinates of cells
# KNN = how many nearest neighbours to find
find_KNN <- function(xy, KNN = 25){

	euc_dist<-as.matrix(dist(xy,upper=T))
	NN<-apply(euc_dist,2,function(x) order(x,decreasing=F)[1:(KNN+1)])
	return(NN)
}

decode_stimulus <- function(log_lik, xy, k, ks){
	n_stim <- dim(log_lik)[2]
	n_reps <- dim(log_lik)[1] / n_stim
	n_cells <- dim(log_lik)[3]
		
	NN2 <- find_KNN(xy, k)
	decoded <- array(dim = c(length(ks), n_reps * n_stim, n_cells))
	for (i in 1 : length(ks)){
		cat(i,"\r")
		for (j in 1 : n_cells){
	
			k_log_lik <- log_lik[, , NN2[1 : ks[i], j]]
			if (ks[i] == 1){
				decoded[i, ,j] <- apply(k_log_lik, 1, which.max)
			}else{
				decoded[i, ,j] <- apply(apply(k_log_lik, c(1, 2), sum), 1, which.max)
			}
		}
	}
    return(decoded)
}

calc_stim_decoding <- function(actual, decoded){
	stim <-	unique(actual)
	perc_decoded_list <- vector("list", length(stim))
	for (i in stim){
		index <- which(actual == i)
		perc_decoded_list[[i]] <- apply(decoded[, index, ], c(1, 3), function(x) sum(x == actual[index]) / length(x))
	}
	return(perc_decoded_list)
}

calc_quad_curve<-function(xy_r){
	x<-xy_r[,1]
	y<-xy_r[,2]

	x2<-x^2
	quad_fit<-lm(y ~ x + x2)
	return(quad_fit)
}

project_xy<-function(xy_r,grid_size){
    quad_fit <- calc_quad_curve(xy_r)
	x_vals<-seq(0,n_pixels,length.out=grid_size)
	y_pred<-predict(quad_fit,list(x=x_vals,x2=x_vals^2))
	xy_curve<-cbind(x_vals,y_pred)
	
	dists<-apply(xy_r,1,function(xy_r) apply(xy_curve,1,function(pnt) dist(rbind(pnt,xy_r))))
	xy_2_pnt<-apply(dists,2,which.min)
}

heatmap_cols=c(rgb(0,1,0),rgb(1,0,1))

add_alpha <- function(alpha=1,cols=heatmap_cols){
	  apply(sapply(cols, col2rgb)/255, 2, function(x)  rgb(x[1], x[2], x[3], alpha=alpha))  
}

#generates a palette of colours for a heatmap from a vector (called internally)
colPalette<-function(vec,alpha_val=1,cols=heatmap_cols,col_depth=100,Viridis=F){
	if(Viridis){
		cols<-viridis(col_depth)[cut(vec,breaks=col_depth,labels=F)]
	}else{
		colPal<-colorRampPalette(colors=add_alpha(alpha_val,cols=cols),alpha=T)
		cols<-colPal(col_depth)[cut(vec,breaks=col_depth,labels=F)]
	}
	return(cols)
}

martincolourscale=c("green4","royalblue4","magenta4","yellow3","turquoise4","orangered4","yellow4","hotpink4","olivedrab4","#f3fa54ff","#1f1c25ff","deeppink4","#176462ff","#54f0faff","#fa6d54ff","#da32daff","#54fac4ff","gray4","#fbf2f9ff","#602646ff","#a96350ff","orange2","#e4eac1ff","#deee82ff","blue4","#203655ff","#989865ff","#f2e7f7ff","#6ef914ff","#f106e3ff","#060cf1ff","#1f1c25ff");


`%inm%` <- function(rows, mat){
	test <- apply(mat, 1, `==`, rows)
  	any(apply(test, 2, all))
}

rotate <- function(x) t(apply(x, 2, rev))

mask_cells <- function(xy, mask){   
	

		x <- floor(xy[, 1])
		x[x == 0] <- 1
		y <- floor(xy[, 2])
		y[y == 0] <- 1

	    x[x > n_pixels] <- 256
	    y[y > n_pixels] <- 256


		cells_vec <- rep(0, n_pixels * n_pixels)
		cells_vec[x + n_pixels * (y - 1)] <- 1
		cells_mat <- matrix(cells_vec, ncol = n_pixels)
		tectal_cells <- mask * cells_mat 
		tectal_cells <- which(tectal_cells == 1, arr.ind = T)
		tectal_cells_index <- apply(cbind(x, y), 1, `%inm%`, tectal_cells)		
        return(tectal_cells_index)

}


