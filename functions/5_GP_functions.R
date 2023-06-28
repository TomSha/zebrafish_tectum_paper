########################################################################
# calculates the joint exponentiated covariance function across two dimensions 
# params:
# Xi = a matrix of xy spatial coordinates (i_pnts x 2D coordinates)
# Xj = a matrix of xy spatial coordinates (j_pnts x 2D coordinates)
# rho1 = a positive real number which defines the hyperparameter for the length scale in the x direction
# rho2 = a positive real number which defines the hyperparameter for the length scale in the y direction
# alpha = a positive real number which defines the hyperparameter for the marginal standard deviation (ie the magnitude of the range of the function) 
# returns:
# cov_mat = an i_pnts x j_pnts covariance matrix 
########################################################################

cov_joint <- function(Xi, Xj, rho1, rho2, alpha){

	cov_mat <- alpha ^ 2 * exp(-1 / (2 * rho1 ^ 2) * outer(Xi[,1], Xj[,1], "-") ^ 2 -1 / (2 * rho2 ^ 2) * outer(Xi[,2], Xj[,2], "-") ^ 2)
	
	return(cov_mat)
}


########################################################################
# calculates the exponentiated covariance function across one dimension
# params:
# Xi = a vector of 1D spatial coordinates, either the x or y coordinate (i_pnts x 1D coordinate)
# Xj = a matrix of 1D spatial coordinates, either the x or y coordinate (j_pnts x 1D coordinate)
# rho = a positive real number which defines the hyperparameter for the length scale
# alpha = a positive real number which defines the hyperparameter for the marginal standard deviation (ie the magnitude of the range of the function) 
# returns:
# cov_mat = an i_pnts x j_pnts covariance matrix 
########################################################################

cov_marg <- function(Xi, Xj, rho, alpha){

	cov_mat <- alpha ^ 2 * exp(-1 / (2 * rho ^ 2) * outer(Xi, Xj, "-") ^ 2 )

	return(cov_mat)
}

########################################################################
# calculates predicted output values at given input values based on input values where the output is known (2 dimensional)
# params:
# X = a matrix of input values (xy spatial coordinates) where the output is known (i_pnts x 2D coordinates)
# X_predict = a matrix of input values (xy spatial coordinates) where you want to predict the output (j_pnts x 2D coordinates)
# Y = a vector of output values which correspond to the X input coordinates (i output values)
# rho1 = a positive real number which defines the hyperparameter for the length scale in the x direction
# rho2 = a positive real number which defines the hyperparameter for the length scale in the y direction
# alpha = a positive real number which defines the hyperparameter for the marginal standard deviation (ie the magnitude of the range of the function) 
# sigma = the noise term in the regression (std. dev)
# returns:
# pred_mean = a vector of the predicted mean of the output at X_predict
# pred_var = a vector of the predicted variance of the output at X_predict
########################################################################

prediction_joint <- function(X, X_predict, Y, rho1, rho2, alpha, sigma){

	n_dat <- nrow(X)
	Sigma_inv <- solve(cov_joint(X, X, rho1, rho2, alpha) + diag(sigma, n_dat))
	Sigma_new <- cov_joint(X_predict, X, rho1, rho2, alpha)
	Sigma_star <- cov_joint(X_predict, X_predict, rho1, rho2, alpha)

	pred_mean <- Sigma_new %*% Sigma_inv %*% Y
	pred_var <- Sigma_star - Sigma_new %*% Sigma_inv %*% t(Sigma_new)
	returnList <- list(pred_mean = pred_mean, pred_var = diag(pred_var))

	return(returnList)
}	

########################################################################
# calculates the predicted output values at given input values based on input values where the output is known (1 dimensional)
# params:
# X = a vector of input values (1D spatial coordinates) where the output is known (i_pnts)
# X_predict = a vector of input values (1D spatial coordinates) where you want to predict the output (j_pnts)
# Y = a vector of output values which correspond to the X input coordinates (i output values)
# rho1 = a positive real number which defines the hyperparameter for the length scale in the x direction
# rho2 = a positive real number which defines the hyperparameter for the length scale in the y direction
# alpha = a positive real number which defines the hyperparameter for the marginal standard deviation (ie the magnitude of the range of the function) 
# sigma = the noise term in the regression (std. dev)
# returns:
# pred_mean = a vector of the predicted mean of the output at X_predict
# pred_var = a vector of the predicted variance of the output at X_predict
########################################################################

prediction_marg <- function(X, X_predict, Y, rho, alpha, sigma){

	n_dat <- length(X)
	Sigma_inv <- solve(cov_marg(X, X, rho, alpha) + diag(sigma, n_dat))
	Sigma_new <- cov_marg(X_predict, X, rho, alpha)
	Sigma_star <- cov_marg(X_predict, X_predict, rho, alpha)

	pred_mean <- Sigma_new %*% Sigma_inv %*% Y
	pred_var <- Sigma_star - Sigma_new %*% Sigma_inv %*% t(Sigma_new)
	returnList <- list(pred_mean = pred_mean, pred_var = diag(pred_var))
	return(returnList)
}

########################################################################
# calculates the predicition_joint() and prediction_marg() for n samples drawn from MCMC inference on the hyperparameters of the exponentiated covariance function
# params:
# samples = a list of vectors where each vector contains n draws from the posterior of a hyperparameter
# X = a matrix of input values (xy spatial coordinates) where the output is known (i_pnts x 2D coordinates)
# X_predict = a matrix of input values (xy spatial coordinates) where you want to predict the output (j_pnts x 2D coordinates)
# Y = a vector of output values which correspond to the X input coordinates (i output values)
# n_samps = the number of samples to use to calculate the regression
#
# returns:
# joint_mean = the predicted mean value of the function over the 2D spatial coordinates in X_predict
# join_var = the predicted variance of the the value of the function over the 2D spatial coordinates in X_predict
# marg1_mean = the predicted mean value of the function over the X spatial dimesion at X_predict[,1]
# marg1_var = the predicted varianace of the value of the function over the X spatial dimention at X_predict[,1]
# marg2_mean = the predicted mean value of the function over the Y spatial dimension at X_predcit[,2]
# marg2_var = the predicted variance  of the value of the function over the Y spatial dimension at X_predict[,2]
########################################################################

calc_reg <- function(samples, X, Y, X_predict, n_samps = 1000){

	n_dat <- nrow(X)

	joint_mean <- matrix(NA, ncol = nrow(X_predict), nrow = n_samps)
	joint_var <- matrix(NA, ncol = nrow(X_predict), nrow = n_samps)

	marg1_mean <- matrix(NA, ncol = nrow(X_predict), nrow = n_samps)
	marg1_var <- matrix(NA, ncol = nrow(X_predict), nrow = n_samps)

	marg2_mean <- matrix(NA, ncol = nrow(X_predict), nrow = n_samps)
	marg2_var <- matrix(NA, ncol = nrow(X_predict), nrow = n_samps)
		

        for(j in 1:n_samps){
			cat(j,"\r")
            i <- n_samps - j+1
			rho1 <- samples$rho1[i]
			rho2 <- samples$rho2[i]
			alpha <- samples$alpha[i]
			sigma <- samples$sigma[i]
			pred_joint <- prediction_joint(X, X_predict, Y, rho1, rho2, alpha, sigma)
			joint_mean[i,] <- pred_joint$pred_mean
			joint_var[i,] <- pred_joint$pred_var

			pred_marg1 <- prediction_marg(X[,1], X_predict[,1], Y, rho1, alpha, sigma)
			marg1_mean[i,] <- pred_marg1$pred_mean
			marg1_var[i,] <- pred_marg1$pred_var

			pred_marg2 <- prediction_marg(X[,2], X_predict[,2], Y, rho2, alpha, sigma)
			marg2_mean[i,] <- pred_marg2$pred_mean
			marg2_var[i,] <- pred_marg2$pred_var
        }

   returnList <- list(joint_mean = joint_mean, joint_var = joint_var, marg1_mean = marg1_mean, marg1_var = marg1_var, marg2_mean = marg2_mean, marg2_var = marg2_var)

   return(returnList)
}

