### distributions.R:  distribution and parameter estimation functions for the lcmix package

### UNIVARIATE DISTRIBUTION DEFINITIONS

## "WeiSD" (Weibull, shape-decay parameterization) functions;
## if x ~ WeiSD(shape, decay) then:
## f(x) = shape * decay * x^(shape-1) * exp(-decay * x^shape)
## F(x) = 1 - exp(-decay * x^shape)
## S(x) = exp(-decay * x^shape)
## H(x) = decay * x^shape
## h(x) = decay * shape * x^(shape-1)
## This implies that the "shape" parameters for WSR and the standard R 
## parameterization of the Weibull distribution are the same, while
## using the "scale" parameter in the standard R parameterization we
## have decay = scale^(-shape), or equivalently, scale = decay^(-1/shape).

# probability density function
dweisd <- function(x, shape, decay, log=FALSE)
	dweibull(x, shape=shape, scale=decay^(-1/shape), log=log)

# cumulative density function
pweisd <- function(q, shape, decay, lower.tail=TRUE, log.p=FALSE)
	pweibull(q, shape=shape, scale=decay^(-1/shape),
			 lower.tail=lower.tail, log.p=log.p)

# quantile function
qweisd <- function(p, shape, decay, lower.tail=TRUE, log.p=FALSE)
	qweibull(p, shape=shape, scale=decay^(-1/shape),
			 lower.tail=lower.tail, log.p=log.p)

# random generation function
rweisd <- function(n, shape, decay)
	rweibull(n, shape=shape, scale=decay^(-1/shape))

# expected value
eweisd <- function(shape, decay) decay^(-1/shape) * gamma(1 + 1/shape)

# variance
vweisd <- function(shape, decay)
	decay^(-2/shape) * gamma(1 + 2/shape) - eweisd(shape, decay)^2

# median
medweisd <- function(shape, decay) decay^(-1/shape) * log(2)^(1/shape)

# for x ~ WeiSD(shape, decay), return the rth raw moment of x, i.e. E(x^r)
rmomweisd <- function(r, shape, decay) decay^(-r/shape) * gamma(1 + r/shape)

### MULTIVARIATE DISTRIBUTION DEFINITIONS

## The categorical distribution.

# density function
dcategorical <- function(x, prob, log=FALSE)
	if(log) log(prob)[x] else prob[x]

# random generation function
rcategorical <- function(n, prob)
	sample(1:length(prob), n, replace=TRUE, prob=prob)

## Multivariate gamma distribution using normal copula (TO DO:  figure out an elegant way to allow for scale parameter; it may be as simple as throwing the "scale=1/rate" default into the function calls and then handing scale off to the dgamma() and pgamma() functions along with shape and rate ... 2011-07-05 DD)

# density function
dmvgamma <- function(X, shape, rate, corr=diag(ncol(X)), log=FALSE)
{	
	## gather statistics

	D = ncol(X)
	
	cdf = sapply(1:D, function(d) pgamma(X[,d], shape[d], rate[d]))
	normcdf = qnorm(cdf)
	normcdf = force.range(normcdf, LC_MINSTDNORM, LC_MAXSTDNORM)
		# deal with outliers so small (large) that their effective CDF is 0 (1); this represents a robustification approach
	
	## calculate copula contribution to log-PDF
	
	res.copula.add = dmvnorm(normcdf, cov=corr, log=TRUE)
	res.copula.sub = rowSums(dnorm(normcdf, log=TRUE))
	res.copula = res.copula.add - res.copula.sub
	
	## calculate marginal contribution to log-PDF
	
	res.data = rowSums(sapply(1:D, function(d)
		dgamma(X[,d], shape[d], rate[d], log=TRUE)))
	
	## final calculations and return
	
	retn = res.copula + res.data
	
	if(log) retn else exp(retn)
}

# random generation function
rmvgamma <- function(N, shape=1, rate=1, corr=diag(length(shape)))
{
	## extract parameters, do sanity checks, deal with univariate case
	
	if(!is.matrix(corr) || !isSymmetric(corr))
		stop("'corr' must be a symmetric matrix")
	D = ncol(corr)
	
	Ds = length(shape)
	if(Ds > D)
		warning("'shape' longer than width of 'corr', truncating to fit")
	if(Ds != D)
		shape = rep(shape, length.out=D)
	
	Dr = length(rate)
	if(Dr > D)
		warning("'rate' longer than width of 'corr', truncating to fit")
	if(Dr != D)
		rate = rep(rate, length.out=D)
	
	if(D == 1) rgamma(N, shape, rate)
	
	## generate standard multivariate normal matrix, convert to CDF
	
	Z = rmvnorm(N, cov=corr)
	cdf = pnorm(Z)
	
	## convert to gamma, return
	
	sapply(1:D, function(d) qgamma(cdf[,d], shape[d], rate[d]))
}

## Multivariate normal distribution

# density function
# TO DO:  can we speed this up (and still retain numerical stability) by using the more familiar LU-decomposition-based expressions such as solve() and determinant() for the terms of the log-PDF?	 Something to experiment with, once we're sure that everything else is working ... 2011-07-05 DD
dmvnorm <- function(X, mean=rep(0, ncol(X)), cov=diag(ncol(X)), log=FALSE)
{
	## preprocessing
	
	# munge data PRN
	if(is.data.frame(X))
		X = as.matrix(X)
	else if(is.vector(X))
		X = matrix(X, nrow=1)
	
	# gather statistics
	N = nrow(X)
	D = ncol(X)
	
	# adjust data by mean
	X = t(t(X) - mean)
	
	## build terms of log-PDF
	
	qr.cov = qr(cov)
	determinant.term = sum(log(abs(diag(qr.cov$qr))))
	
	constant.term = D * LC_LOG2PI
	
	distance.term = if(N == 1)
		as.vector(X %*% qr.solve(qr.cov, t(X)))
	else
		rowSums((X %*% qr.solve(qr.cov)) * X)
	
	## calculate and return (log-)density
	
	res = -0.5 * (determinant.term + constant.term + distance.term)
	
	if(log) res else exp(res)
}

# random generation function
rmvnorm <- function(N, mean=NULL, cov=NULL) 
{
	## munge parameters PRN and deal with the simplest univariate case
	
	if(is.null(mean))
		if(is.null(cov))
			return(rnorm(N))
		else
			mean = rep(0, nrow(cov))
	else if (is.null(cov))
		cov = diag(length(mean))
	
	## gather statistics, do sanity checks
	
	D = length(mean)
	if (D != nrow(cov) || D != ncol(cov)) 
		stop("length of mean must equal nrow and ncol of cov")
	
	E = eigen(cov, symmetric=TRUE)
	if (any(E$val < 0)) 
		stop("Numerically negative definite covariance matrix")
	
	## generate values and return
	
	mean.term = mean
	covariance.term = E$vec %*% (t(E$vec) * sqrt(E$val))
	independent.term = matrix(rnorm(N*D), nrow=D)
	
	t(mean.term + (covariance.term %*% independent.term))
}

## Multivariate Weibull (WeiSD) distribution using normal copula

# density function (X is a N-by-D matrix, shape and decay are D-length vectors, corr is a D-by-D matrix)
dmvweisd <- function(X, shape, decay, corr=diag(ncol(X)), log=FALSE)
{	
	## gather statistics

	D = ncol(X)
	
	cdf = sapply(1:D, function(d) pweisd(X[,d], shape[d], decay[d]))
	normcdf = qnorm(cdf)
	normcdf = force.range(normcdf, LC_MINSTDNORM, LC_MAXSTDNORM)
		# deal with outliers so small (large) that their effective CDF is 0 (1); this represents a robustification approach
	
	## calculate copula contribution to log-PDF
	
	res.copula.add = dmvnorm(normcdf, cov=corr, log=TRUE)
	res.copula.sub = rowSums(dnorm(normcdf, log=TRUE))
	res.copula = res.copula.add - res.copula.sub
	
	## calculate marginal contribution to log-PDF
	
	res.data = rowSums(sapply(1:D, function(d)
		dweisd(X[,d], shape[d], decay[d], log=TRUE)))
	
	## final calculations and return
	
	retn = res.copula + res.data
	
	if(log) retn else exp(retn)
}

# random generation function
rmvweisd <- function(N, shape=1, decay=1, corr=diag(length(shape)))
{
	## extract parameters, do sanity checks, deal with univariate case
	
	if(!is.matrix(corr) || !isSymmetric(corr))
		stop("'corr' must be a symmetric matrix")
	D = ncol(corr)
	
	Ds = length(shape)
	if(Ds > D)
		warning("'shape' longer than width of 'corr', truncating to fit")
	if(Ds != D)
		shape = rep(shape, length.out=D)
	
	Dd = length(decay)
	if(Dd > D)
		warning("'decay' longer than width of 'corr', truncating to fit")
	if(Dd != D)
		decay = rep(decay, length.out=D)
	
	if(D == 1) rweisd(N, shape, decay)
	
	## generate standard multivariate normal matrix, convert to CDF
	
	Z = rmvnorm(N, cov=corr)
	cdf = pnorm(Z)
	
	## convert to Weibull (WeiSD), return
	
	sapply(1:D, function(d) qweisd(cdf[,d], shape[d], decay[d]))
}

### PARAMETER ESTIMATION FUNCTIONS, WITH OPTIONAL WEIGHTS

## gamma distribution

thetahat.gamma <- function(x, w=1, aslist=TRUE)
{
	## sanity check and gather statistics
	
	mx = min(x)
	if(mx < 0)
		stop("data must be non-negative")
	if(mx < LC_EPS)
		x = force.range(x, LC_EPS) # prevent instability
	
	N = length(x)
	Nw = length(w)
	if(Nw > N)
		warning("weights longer than data, truncating to fit")
	if(Nw != N)
		w = rep(w, length.out=N)
	
	wsum = sum(w)
	lx = log(x)
	
	## take advantage of existing multivariate normal estimators to do initial method-of-moments parameter estimation
	
	xmoments = thetahat.norm(x, w)
	shape.guess = xmoments$mean^2 / xmoments$var
	
	## profile likelihood parameter estimation
	
	toptim <- function(shape)
	{
		rate = shape*wsum / sum(w*x)
		
		term1 = shape * log(rate) * wsum
		term2 = lgamma(shape) * wsum
		term3 = (shape-1) * sum(w * lx)
		term4 = shape * wsum # = rate * sum(w*x)
		
		-(term1 - term2 + term3 - term4)
	}
	
	shape = nlminb(shape.guess, toptim, lower=LC_EPS)$par
	rate = shape*wsum / sum(w*x)
	scale = 1/rate
	
	## package it up, send it back
	
	if(aslist)
		list(shape=shape, rate=rate, scale=scale)
	else
		c(shape=shape, rate=rate, scale=scale)
}

thetahat.mvgamma <- function(X, w=1, aslist=TRUE)
{
	## munge data PRN, gather statistics, do sanity checks
	
	if(!is.matrix(X)) X = as.matrix(X)
	
	N = nrow(X)
	D = ncol(X)
	Nw = length(w)
	if(Nw > N)
		warning("weights longer than data, truncating to fit")
	if(Nw != N)
		w = rep(w, length.out=N)
	
	## marginal parameter estimation
	
	res = sapply(1:D, function(d) thetahat.gamma(X[,d], w, FALSE))
	colnames(res) = colnames(X)
	
	retn = list(shape=res["shape",], rate=res["rate",])

	## copula parameter estimation
	
	cdf = sapply(1:D, function(d) pgamma(X[,d], retn$shape[d], retn$rate[d]))
	normcdf = qnorm(cdf)
	normcdf = force.range(normcdf, LC_MINSTDNORM, LC_MAXSTDNORM)
		# deal with outliers so small (large) that their effective CDF is 0 (1); this represents a robustification approach
	
	retn$corr = cov2cor(t(normcdf) %*% (w*normcdf) / sum(w))
		# Phi-transformed CDFs are already assumed to have mean 0, so no need to estimate means and center; "cov2cor" is a hack to make sure the result is a valid correlation matrix, and fixing it is TO DO (2011-07-16 DD)
	
	## package it up, send it back
	
	if(aslist)
		return(retn)
	else
		unlist(retn)
}

## normal distribution

thetahat.norm <- function(x, w=1, aslist=TRUE)
{
	## gather statistics, do sanity checks
	
	N = length(x)
	Nw = length(w)
	if(Nw > N)
		warning("weights longer than data, truncating to fit")
	if(Nw != N)
		w = rep(w, length.out=N)
	
	## parameter estimation
	
	wsum = sum(w)
	
	mean = sum(w * x) / wsum
	var = sum(w * x^2) / wsum - mean^2
	sd = sqrt(var)
	
	## package it up, send it back
	
	if(aslist)
		list(mean=mean, var=var, sd=sd)
	else
		c(mean=mean, var=var, sd=sd)
}

thetahat.mvnorm <- function(X, w=1, aslist=TRUE)
{
	## munge data PRN, gather statistics, do sanity checks
	
	if(!is.matrix(X)) X = as.matrix(X)
	
	N = nrow(X)
	Nw = length(w)
	if(Nw > N)
		warning("weights longer than data, truncating to fit")
	if(Nw != N)
		w = rep(w, length.out=N)
	
	## parameter estimation
	
	wsum = sum(w)
	wX = w * X
	
	mean = colSums(wX) / wsum
	cov = t(X) %*% (wX) / wsum - tcrossprod(mean)
		# "tcrossprod(mean)" is equivalent to "mean %*% t(mean)" or "outer(mean, mean)" but should be slightly faster
	
	## package it up, send it back
	
	if(aslist)
		list(mean=mean, cov=cov)
	else
		rbind(mean, cov)
}

## Weibull (WeiSD) distribution

thetahat.weisd <- function(x, w=1, aslist=TRUE)
{	
	## sanity check and gather statistics
	
	mx = min(x)
	if(mx < 0)
		stop("data must be non-negative")
	if(mx < LC_EPS)
		x = force.range(x, LC_EPS) # prevent instability
	
	N = length(x)
	Nw = length(w)
	if(Nw > N)
		warning("weights longer than data, truncating to fit")
	if(Nw != N)
		w = rep(w, length.out=N)
	
	lx = log(x)
	wsum = sum(w)
	swlx = sum(w * lx)
	
	## initial method-of-moments parameter estimation
	
	shape.guess = pi / sqrt(6 * thetahat.norm(lx, w)$var)
		# see Lawless (1982) pp. 18-19
	
	## profile likelihood parameter estimation
	
	toptim <- function(shape)
		-log(shape)*wsum + log(sum(w * x^shape))*wsum - (shape-1)*swlx
	
	shape = nlminb(shape.guess, toptim, lower=LC_EPS)$par
	decay = wsum / sum(w * x^shape)
	
	## package it up, send it back
	
	if(aslist)
		list(shape=shape, decay=decay)
	else
		c(shape=shape, decay=decay)
}

thetahat.mvweisd <- function(X, w=1, aslist=TRUE)
{
	## munge data PRN, gather statistics, do sanity checks
	
	if(!is.matrix(X)) X = as.matrix(X)
	
	N = nrow(X)
	D = ncol(X)
	Nw = length(w)
	if(Nw > N)
		warning("weights longer than data, truncating to fit")
	if(Nw != N)
		w = rep(w, length.out=N)
	
	## marginal parameter estimation
	
	res = sapply(1:D, function(d) thetahat.weisd(X[,d], w, FALSE))
	colnames(res) = colnames(X)
	
	retn = list(shape=res["shape",], decay=res["decay",])

	## copula parameter estimation
	
	cdf = sapply(1:D, function(d) pweisd(X[,d], retn$shape[d], retn$decay[d]))
	normcdf = qnorm(cdf)
	normcdf = force.range(normcdf, LC_MINSTDNORM, LC_MAXSTDNORM)
		# deal with outliers so small (large) that their effective CDF is 0 (1); this represents a robustification approach
	
	retn$corr = cov2cor(t(normcdf) %*% (w*normcdf) / sum(w))
		# Phi-transformed CDFs are already assumed to have mean 0, so no need to estimate means and center; "cov2cor" is a hack to make sure the result is a valid correlation matrix, and fixing it is TO DO (2011-07-16 DD)
	
	## package it up, send it back
	
	if(aslist)
		return(retn)
	else
		unlist(retn)
}
