### mixmod.R:  define the mixmod (single data mixture models) function for the lcmix package.

mixmod <- function(X, K, family=names(LC_FAMILY), iter.max=LC_ITER_MAX)
{
	## munge data, gather preliminary information, and do sanity checks
	
	dname = deparse(substitute(X))
	dattr = attributes(X)
	
	X = drop(as.matrix(X)) # get matrix or vector form
	if(is.matrix(X)) {
		N = nrow(X)
		D = ncol(X)
		univariate = FALSE
	} else {
		N = length(X)
		D = 1
		univariate = TRUE
	}
	
	family = match.arg(family)
	dummy = .sanity_check_data(X, family) # if we get through this, we're okay
	distn = ifelse(univariate, LC_FAMILY[[family]]$uni, LC_FAMILY[[family]]$multi)
	
	kvec = 1:K
	
	npar.hidden = K-1
	npar.observed = K * .npar_observed(X, distn)
	npar = npar.hidden + npar.observed
	bicsub = npar * log(N)
	
	## define internal functions
	
	# weight matrices given PDFs and previous iteration's parameters
	estep <- function(pdfs, params)
	{
		numerator = t(t(pdfs$G) * params$hidden$prob)
		W = as.matrix(numerator / pdfs$fX)
			# "as.matrix" to deal with the pathological K=1 case
		
		list(W=W)
	}
	
	# parameter estimates given weight matrices
	mstep.observed <- .mstep_observed(distn)
	mstep <- function(weights)
	{
		list(hidden = list(prob=colMeans(weights$W)),
		     observed = mstep.observed(X, kvec, weights$W))
	}
	
	# PDFs given parameter estimates
	pstep.observed <- .pstep_observed(distn)
	pstep <- function(params)
	{	
		G = pmax(pstep.observed(X, kvec, params$observed), LC_EPSEPS)
			# "pmax" to deal with numerical instability of outliers; we use LC_EPSEPS here because LC_EPS is way too large for this purpose (and LC_EPSEPS is still large enough to avoid triggering badness)
		
		fX = drop(G %*% params$hidden$prob) # = f(X)
			# "drop" so we have an N-length vector rather than a N-by-1 matrix
		
		list(G=G, fX=fX)
	}
	
	# iteration statistics given iteration and results of previous steps
	sstep <- function(iter, weights, params, pdfs)
	{
		llik = sum(log(pdfs$fX))
			# log-likelihood
		qval = sum(.qfun_mixmod(weights, params, pdfs))
		bic = 2*llik - bicsub
		
		c(iter=iter, llik=llik, qval=qval, bic=bic)
	}
	
	## initialize
	
	# quasi E-step
	if(family %in% LC_NONNEGFAM) # convert to (-Inf, Inf) for Gaussian weights
		weights = list(W = .init_weights(log(X), K))
	else
		weights = list(W = .init_weights(X, K))
		
	# M-step
	params = mstep(weights)
		
	# post-EM calculations
	pdfs = pstep(params)
	stats = sstep(0, weights, params, pdfs)
	
	# save results
	iteration.params = list(params)
	iteration.stats = list(stats)
	
	## iterate
	
	old_llik = stats[["llik"]] - 1 # dummy value to ensure iteration
	
	iter = 1
	
	while(abs(stats[["llik"]]/old_llik - 1) > LC_EPS && iter <= iter.max)
	{
		# preserve previous log-likelihood
		old_llik = stats[["llik"]]
		
		# E-step
		weights = estep(pdfs, params)
		
		# M-step
		params = mstep(weights)
		
		# post-EM calculations
		pdfs = pstep(params)
		stats = sstep(iter, weights, params, pdfs)
		
		# update iteration count and save results
		iter = iter + 1
		iteration.params[[iter]] = params
		iteration.stats[[iter]] = stats
			# increment iter before saving because iteration.params[[1]] is actually the parameters from iteration 0, etc., and the same for iteration.params	
	}
	
	## final calculations
	
	weights = estep(pdfs, params)
	
	posterior = weights$W
	
	assignment = apply(posterior, 1, which.max)
	
	## package and return
	
	iteration.stats = Reduce(rbind, iteration.stats)
	rownames(iteration.stats) = NULL
	iteration.stats = as.data.frame(iteration.stats)
	
	retn = list(N                = N,
	            D                = D,
	            K                = K,
	            X                = X,
	            npar             = npar,
	            npar.hidden      = npar.hidden,
	            npar.observed    = npar.observed,
	            iter             = iter,
	            params           = params,
	            stats            = stats,
	            weights          = weights,
	            pdfs             = pdfs,
	            posterior        = posterior,
	            assignment       = assignment,
	            iteration.params = iteration.params,
	            iteration.stats  = iteration.stats,
	            family           = family,
	            distn            = distn,
	            iter.max         = iter.max,
	            dname            = dname,
	            dattr            = dattr)
	class(retn) = "mixmod"
	
	return(retn)
}
