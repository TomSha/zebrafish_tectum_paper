
PlotPrior <- function(mu0,tau0,alpha0,beta0,mur,tr,data=NULL){
		grid=expand.grid(mur,tr)
		f=apply(grid,1,function(x) x[2]^(alpha0-0.5)*exp(-x[2]/2*(tau0*(x[1]-mu0)^2+2*beta0))*beta0^alpha0/gamma(alpha0)*sqrt(tau0/2/pi))
		f=f/sum(f)
		layout(matrix(1:2,1))
		image(x=mur,y=tr,z=matrix(f,length(mur)),col=gray.colors(1000))
		if(!is.null(data)) points(data[,1],data[,2],col="red")
		contour(x=mur,y=tr,z=matrix(f,length(mur)),add=T)
		#plot(tr,tr^(alpha0-0.5)*exp(-tr/2*(tau0*(mu0)^2+2*beta0))*beta0^alpha0/gamma(alpha0)*sqrt(tau0/2/pi),type='l',ylab="P(t|m=0)")
		plot(tr,dgamma(tr,alpha0,beta0),type='l',ylab="P(t)")

}

colMax <- function(X) apply(X,2,max)
rowMax <- function(X) apply(X,1,max)

########################################################################
# calculates the marginal log likelihood of a gaussian with unknown mean and std. dev
# params:
# data = a numeric vector  
# mu0 = prior mean
# tau0 = prior precision
# alpha0, beta0 = parameters of the gamma prior on the precision
# This set of parameters is good
# mu0=0,tau0=.01,alpha0=.5,beta0=10
# returns:
# logML = log marginal likelihood of the data
########################################################################

MarginalLogLikelihood <- function(data,mu0=0,tau0=0.01,alpha0=0.5,beta0=1){
		
	em=mean(data)
	n=length(data)
	tau_n=tau0+n
	alpha_n=alpha0+n/2
	beta_n=beta0+0.5*sum((data-em)^2)+0.5*tau0/tau_n*n*(em-mu0)^2

        logML=lgamma(alpha_n)-lgamma(alpha0)+alpha0*log(beta0)-alpha_n*log(beta_n)+0.5*log(tau0/tau_n)-n/2*log(2*pi)
	 
	return(logML)
}



########################################################################
# will calculate the marg loglikelihood for every possible partition of n stimuli - a partition of a set n is defined as a set of nonempty, pairwise disjoint subsets of n whose union is n.
# params:
# dat = vector of responses of one cell, ordered by stimulus eg {stim1_rep1, stim1_rep2, stim2_rep1, stim2_rep2,...}
# n_stim = number of stimuli
# n_rep = number of repetitions of each stimulus
# parti = output from listParts() which gives a list of every possible partition of a set of size n
# returns:
# modelML = vector of the loglikelihoods for each partition of the data
########################################################################

calc_loglikelihood<-function(dat, n_stim, n_rep, parti){

	n_parts<-length(parti)
	modelML<-rep(NA, n_parts)
	indicator<-rep(1 : n_stim, each = n_rep)


	for(i in 1 : n_parts){
		cat(i, "of", n_parts, "\r")
		k <- length(parti[[i]])
		lvec <- rep(NA,k)
		for(j in 1 : k){
			    lvec[j] <- MarginalLogLikelihood(dat[indicator %in% parti[[i]][[j]]])
		}
		modelML[[i]] <- sum(lvec)
	}
	return(modelML)
}

########################################################################
# calculates the response selectivity of a cell for cells whose responses are best described by partitioning them into 2 subsets (as determined by the loglikelihood). If the responses are not best described by 2 sets the function returns NA. The selectivity is bounded between -1 and +1.
# params:
# n = model number which best describes the cell
# r = response of a cell (2D array [mean from each repetition,stim])
# parti = list of partitions of the stimuli
# returns:
# selec = the response selectivity
########################################################################

calc_model_selectivity <- function(n, r, parti){
	if(length(parti[[n]])==2){
		gauss1 <- parti[[n]]$`2`
		gauss2 <- parti[[n]]$`1`
		
		mu1 <- mean(r[,gauss1])
		mu2 <- mean(r[,gauss2])

		mu1[mu1<0] <- 0
		mu2[mu2<0] <- 0

		selec <- (mu1 - mu2) / (mu1 + mu2)
		return(selec)
	}else{
		return(NA)
	}
}
	
########################################################################
# Assigns neurons to their subtype
# 
# params:
# responses = a matrix of neuronal responses to stimuli (neurons x response)
# epochs = a vector of stimuli labels
#
# returns:
# model_n = a vector of which subtype a neuron has been associated to
# parti = a list of the partitions of the stimuli into non-empty subsets
# LL = a matrix of the marginal loglikelihood for each possible subtype (subtype x neurons)
# model_2gauss = a vector of whether the subtype partitions the data into 2 subsets
########################################################################

cluster_cells <- function(responses, epochs){

	n_stim <- length(epochs)
	n_rep <- ncol(responses)/n_stim

	# generate a list of all possible partitions of the set of stimuli
	parti <- listParts(n_stim)
	n_parts <- length(parti)

	# for every cell and each patition calculate the likelihood of a cell's response
	LL <- apply(responses, 1, calc_loglikelihood, n_stim, n_rep, parti)

	# for each cell calculate which partition (model) has the maximum likelihood
	model_n <- apply(LL, 2, which.max)

	# find which models split the stimuli into two sets eg {stim A} {stim B..D} or {stim A, stim D} {stim B stim C} etc
	model_length <- sapply(parti, length)
	model_2gauss <- which(model_length == 2)

	returnList <- list(model_n = model_n, parti = parti, LL = LL, model_2gauss = model_2gauss)
		
}
	

