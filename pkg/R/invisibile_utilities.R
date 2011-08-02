### invisible_utilities:  invisible utility functions for the lcmix package

# Calculate the initial weight matrix ("W") for data "X" with "K" means, means being chosen appropriately. The matrix that is returned is N-by-K.
.init_weights <- function(X, K)
{
	X = scale(as.matrix(X))
	
	if(K == 1) # pathological case
		matrix(1, nrow=nrow(X), ncol=1)
	else
		.init_weights_from_means(X, uniform.quantiles(X, K, decreasing=TRUE))
}

# Calculate the initial weight matrix ("W") from a set of K means, "mu", which must be of the same width (number of columns) as "X", i.e. X is a vector of length N and mu is a vector of length K, or X is a N-by-D matrix and mu is a K-by-D matrix, where K < N.  The matrix that is returned is N-by-K.
.init_weights_from_means <- function(X, mu)
{
	X = scale(as.matrix(X))
	mu = as.matrix(mu)
	K = nrow(mu)
	
	Dsq = sapply(1:K, function(k) euclidean(X, mu[k,]))
	
	iDsq = 1 / (Dsq + LC_EPS)
	retn = iDsq / rowSums(iDsq)
	
	mu = (t(retn) %*% X) / colSums(retn)
	ranks = apply(-mu, 2, rank)
		# "-mu" because we want greatest to least
	rowranks = apply(ranks, 1, prod)
	retn = retn[,order(rowranks)]
	
	return(retn)
}

# Workhorse function for (unweighted) MLE of covariance; takes a D-by-N matrix, a D-length vector of means, and N, and returns a MLE covariance matrix.
.mleTcov <- function(tX, mu, N) { tXadj = tX - mu ; tcrossprod(tXadj) / N }

# Return a function which performs weighted parameter estimation for the single-data mixture model M-step (can also be used in the multiple-data mixture model so long as no parameters are shared, if handled appropriately.)
.mstep_observed <- function(distn)
{
	switch(distn,
	
		norm = function(X, kvec, W)
		{
			res = sapply(kvec, function(k)
				thetahat.norm(X, W[,k], aslist=FALSE))
			list(mean=res["mean",], var=res["var",], sd=res["sd",])
				# the kth element of $mean is the mean for the kth component, and similarly for $var and $sd
		},
		
		mvnorm = function(X, kvec, W)
		{
			res = sapply(kvec, function(k)
				thetahat.mvnorm(X, W[,k], aslist=TRUE), simplify=FALSE)
			mu = t(sapply(kvec, function(k) res[[k]]$mean))
				# the kth row of mu is the vector of the column means of X for the kth component
			Sigma = lapply(kvec, function(k) res[[k]]$cov)
				# the kth list element of Sigma is the covariance matrix for the kth component
			list(mean=mu, cov=Sigma)
		},
		
		weisd = function(X, kvec, W)
		{
			res = sapply(kvec, function(k)
				thetahat.weisd(X, W[,k], aslist=FALSE))
			list(shape=res["shape",], decay=res["decay",])
				# the kth element of $shape is the shape for the kth component, and similarly for $decay
		},
		
		mvweisd = function(X, kvec, W)
		{
			res = sapply(kvec, function(k)
				thetahat.mvweisd(X, W[,k], aslist=TRUE), simplify=FALSE)
			sigma = t(sapply(kvec, function(k) res[[k]]$shape))
				# the kth row of sigma is the vector of the column shapes of X for the kth component
			delta = t(sapply(kvec, function(k) res[[k]]$decay))
				# the kth row of delta is the vector of the column decays of X for the kth component
			rho = lapply(kvec, function(k) res[[k]]$corr)
			list(shape=sigma, decay=delta, corr=rho)
		},
		
		gamma = function(X, kvec, W)
		{
			res = sapply(kvec, function(k)
				thetahat.gamma(X, W[,k], aslist=FALSE))
			list(shape=res["shape",], rate=res["rate",])
				# the kth element of $shape is the shape for the kth component, and similarly for $rate
		},
		
		mvgamma = function(X, kvec, W)
		{
			res = sapply(kvec, function(k)
				thetahat.mvgamma(X, W[,k], aslist=TRUE), simplify=FALSE)
			sigma = t(sapply(kvec, function(k) res[[k]]$shape))
				# the kth row of sigma is the vector of the column shapes of X for the kth component
			lambda = t(sapply(kvec, function(k) res[[k]]$rate))
				# the kth row of lambda is the vector of the column rates of X for the kth component
			rho = lapply(kvec, function(k) res[[k]]$corr)
			list(shape=sigma, rate=lambda, corr=rho)
		}
	)
}

# Calculate the number of parameters for distribution of observed data in a single-data mixture model
.npar_observed <- function(X, distn)
{
	if(is.matrix(X)) D = ncol(X)
	
	switch(distn,
		norm    = 2,
		mvnorm  = D + D*(D+1)/2,
		weisd   = 2,
		mvweisd = 2*D + D*(D-1)/2,
		gamma   = 2,
		mvgamma = 2*D + D*(D-1)/2)
}

# Return a nicely formatted string version of a family name, suitable for printing
.pretty_family <- function(family)
{
	if(family == "weibull") "Weibull"
	else if(family == "iweibull") "ind. Weibull"
	else if(family == "igamma") "ind. gamma"
	else if(family == "inorm") "ind. normal"
	else family
}

# Return a function which performs PDF estimation for the observed data given the components of a single-data mixture model (can also be used in multiple-data mixture models if handled appropriately.)
.pstep_observed <- function(distn)
{
	switch(distn,
		
		norm = function(X, kvec, obspar) 
			sapply(kvec, function(k)
				dnorm(X, obspar$mean[k], obspar$sd[k])),
		mvnorm = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvnorm(X, obspar$mean[k,], obspar$cov[[k]])),
		
		weisd = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dweisd(X, obspar$shape[k], obspar$decay[k])),
		mvweisd = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvweisd(X, obspar$shape[k,], obspar$decay[k,],
					obspar$corr[[k]])),
		
		gamma = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dgamma(X, obspar$shape[k], obspar$rate[k])),
		mvgamma = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvgamma(X, obspar$shape[k,], obspar$rate[k,],
					obspar$corr[[k]]))
	)
}

# "Q-function" (expected log-likelihood) for single-data (.qfun_mixmod()) and multiple-data (.qfun_mdmixmod()) mixture models.
.qfun_mixmod <- function(weights, params, pdfs)
{
	c(hidden=sum(weights$W %*% log(params$hidden$p)),
	  observed=sum(weights$W * log(pdfs$G)))
}
.qfun_mdmixmod <- function(weights, params, pdfs)
{
	zvec = names(weights$V)
	
	term1 = sum(weights$U %*% log(params$hidden$prob0))
	term2 = sum(sapply(zvec, function(z)
				sum(sapply(weights$V[[z]], function(V_z_elem)
					sum(V_z_elem %*% t(log(params$hidden$cprob[[z]])))))))
	term3 = sum(sapply(zvec, function(z)
				sum(weights$W[[z]] * log(pdfs$G[[z]]))))
	
	c(hidden=term1+term2, observed=term3)
}
		
# Ensure that input data to a single-data mixture model (or a list component in a multiple-data mixture model) meets the requirements of the given distribution family.
.sanity_check_data <- function(X, family)
{
	if((family %in% LC_NONNEGFAM) && min(X) <= 0)
	{
		family = .pretty_family(family)
		stop(sprintf("%s family requires non-negative data", family))
	}
	
	return(0) # normal result, got through the above
}