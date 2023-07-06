
########################################################################
# Code from Giovanni Diana. See Diana et al (2019) Bayesian inference of neuronal assemblies for details.
# Decomposes fluorescent timeseries data into the sum of calcium transient (ck), baseline activity (bk) and a source of Gaussian noise. 
#
# params:
# y = vector of fluorescent timeseries data
# wdt = bernoulli rate
# lam = probe decay
# varB = variance of the baseline
# varC = variance of the signal
# B0 = initial baseline value
# Cmean = initial calcium transient activity value
# freq = est of spike frequency
#
# returns:
# c = vector of calcium signal timeseries
# b = vector of baseline signal timeseries
# sks = vector of binarised output of activity 
# loglik = log-likelihood 
########################################################################

bcl<- function(y,wdt,lam,varB,varC,B0,Cmean,freq=5,mode="GAUS"){
	varX=var(diff(y)[diff(y)<0]);
	N=length(y);
	B=rep(NA,N);
	c=rep(0,N);
	sks=rep(0,N);
    B[1]=B0;
	loglik=0;
	dt=1./freq;
	for(t in 2:N){
		cnew=c[t-1]*exp(-lam*dt);
		Bnew=(varX*B[t-1]+varB*dt*(y[t]-cnew))/(varX+varB*dt);
		if(mode=="GAUS"){
			logp0=log(1-wdt)-0.5*log(2*pi)-0.5*log(varX+varB*dt)-(y[t]-cnew-B[t-1])^2/(2*varX+2*varB*dt);
			cspike=Cmean+cnew+(y[t]-cnew-B[t-1])/(1+varB*dt/varC+varX/varC)
			Bspike=B[t-1]+varB*dt/varC*(cspike-cnew-Cmean);
			logp1=log(wdt)-0.5*log(2*pi)-0.5*log(varX+varB*dt+varC)-(y[t]-cnew-B[t-1]-Cmean)^2/(2*varX+2*varB*dt+2*varC);
		} else if (mode=="EXP"){
			logp0=log(1-wdt)-0.5*log(2*pi)-0.5*log(varX+varB*dt)-(y[t]-cnew-B[t-1])^2/(2*varX+2*varB*dt);
			cspike=y[t]-B[t-1]-varC*(varB*dt+varX);
			Bspike=B[t-1]+varC*varB*dt;
			logp1=log(wdt)-0.5*log((2*pi)^2)-1/2*log(varX)-1/2*log(varB*dt)+log(varC)-(Bspike-B[t-1])^2/(2*varB*dt)*(1+varX/varB/dt) - varC*(cspike-cnew) ;

		} else if (mode=="EXP2"){
			prevs=sks[t-1]+1;
			logp0=log(1-wdt[prevs])-0.5*log(2*pi)-0.5*log(varX+varB*dt)-(y[t]-cnew-B[t-1])^2/(2*varX+2*varB*dt);
			cspike=y[t]-B[t-1]-varC*(varB*dt+varX);
			Bspike=B[t-1]+varC*varB*dt;
			logp1=log(wdt[prevs])-0.5*log((2*pi)^2)-1/2*log(varX)-1/2*log(varB*dt)+log(varC)-(Bspike-B[t-1])^2/(2*varB*dt)*(1+varX/varB/dt) - varC*(cspike-cnew) ;
		}
		
		if(mode=="GAUS"){
			if(logp1<logp0){
				c[t]=cnew;
				B[t]=Bnew;
				loglik=loglik+logp0
			} else {
				c[t]=cspike;
				B[t]=Bspike;
				sks[t]=1;
				loglik=loglik+logp1
			}
		} else if(mode=="EXP" || mode=="EXP2"){
			if(logp1<logp0 || cspike-cnew<=0){
				c[t]=cnew;
				B[t]=Bnew;
				loglik=loglik+logp0;
			} else {
				c[t]=cspike;
				B[t]=Bspike;
				sks[t]=1;
				loglik=loglik+logp1
			}
		}
	}
	
	return(list(c,B,sks,loglik));
}


########################################################################
# Calculates the correlation coefficient between the timepoints of when a 
# neuron was firing and when a visual stimulus was being shown.
# This empirical correlation coefficient is then compared with the correlation 
# coefficient assuming the neuron was firing randomly.
#
# params:
# timelog: an n x 2 data.frame where the 1st column is the time in seconds of when a stimulus started and finished and the 2nd column is the name of the stimulus
# spikes: a n_cells x n_frames matrix of 0s and 1s where 0 indicates a neuron was silent and 1 indicates a neuron was firing
# delay: the delay between starting the imaging session and starting showing the stimulus (usually 0)
# ac_freq: the acquisition frequency of the imaging session
# thresh_val: how many std devs above the random model must a neuron be to be classed as visually responsive
#
# returns:
# a list $corvec is a vector of each neuron's empiricial correlation coefficient 
#	 $threshMean is a vector of the mean correlation coefficients assuming each neuron was firing randomly
#	 $threshSD is a vector of the standard deviations of the correlation coefficient assuming random firing
#	 $thresh_index is a logical vector of whether each neuron is classed as being visually responsive 
# 	 i.e the firing is correlated with the presentation of the visual stimuli
########################################################################

calculate_correlation <- function(timelog, spikes, delay, ac_freq, thresh_val = 0.02){

	# convert timelog from seconds to frames
	timelog[, 1] <- floor((timelog[, 1] + delay) * ac_freq)

	nEP <- nrow(timelog) / 2
	epochs <- timelog[, 2][1:nEP*2]	
	nz = ncol(spikes);
	ncells = nrow(spikes);
	ep <- timelog[, 1];
	matep <- matrix(ep, ncol = 2, byrow = T);

	# adding in 5 seconds after stimulus ends to the stimulus block to allow for the decay of the probe
	matep[, 2] <- matep[, 2] + ceiling(5 * ac_freq)
	corvec <- rep(0, ncells)
	threshMean <- rep(0, ncells)
	threshSD <- rep(0, ncells)
	stim <- rep(-1, nz)
	for(response in 1:nrow(matep)) stim[matep[response, 1] : matep[response, 2]] = 1
	
	# set timepoints to when the neurons aren't firing to -1
	spikes[spikes == 0] <- (-1)
	stim_mat <- matrix(rep(stim, ncells), nrow = ncells, byrow = T)

	# calculate the empirical correlation
	corvec <- rowMeans(spikes * stim_mat)

	# calculate the correlation assuming the neuron fires randomly
	threshMean <- rowMeans(spikes) * rowMeans(stim_mat)
	threshSD <- apply(spikes, 1, sd)
	thresh <- (corvec > (threshMean + (thresh_val * threshSD)))

	returnList <- list(corvec = corvec, threshMean = threshMean, threshSD = threshSD, thresh = thresh)
	return(returnList)
	
}
