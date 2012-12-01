### mixmod.R:  define the mixmod (single data mixture models) function for the lcmix package.

mixmod <- function(X, K, family=names(LC_FAMILY), prior=NULL, 
	iter.max=LC_ITER_MAX, dname=deparse(substitute(X)))
{
	## munge data, gather preliminary information, and do sanity checks
	
	dname = dname
		# keep the implicit call to deparse(substitute(X)) from getting screwed up by whatever manipulation we do below (this is a nasty side-effect of lazy evaluation, I suspect)
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
	distn = ifelse(univariate, LC_FAMILY[[family]]$uni,
		LC_FAMILY[[family]]$multi)
		
	kvec = 1:K
	
	npar.hidden = K-1
	npar.observed = .npar_observed(distn, K, D)
	npar = npar.hidden + npar.observed
	bicsub = npar * log(N)
	
	## define internal functions
	
	estep <- .estep_mixmod()
	mstep <- .mstep_mixmod(distn)
	pstep <- .pstep_mixmod(distn)
	sstep <- .sstep_mixmod(bicsub)
	
	## initialize
	
	# quasi E-step
	weights = list(W = .init_weights(X, K, prior=prior))
	
	# M-step
	params = mstep(X, weights, prior=prior)
				
	# post-EM calculations
	pdfs = pstep(X, params, kvec)
	stats = sstep(0, weights, params, pdfs) # MARK MARK
	
	# save results
	iteration.params = list(params)
	iteration.stats = list(stats)
	
	## iterate
	
	old_llik = stats[["llik"]] - 1 # dummy value to ensure iteration
	
	iter = 1
	
	while(abs(stats[["llik"]]/old_llik - 1) > LC_ITER_TOL && iter <= iter.max)
	{
		# preserve previous log-likelihood
		old_llik = stats[["llik"]]
		
		# E-step
		weights = estep(pdfs, params)
		
		# M-step
		params = mstep(X, weights, params, prior=prior)
		
		# post-EM calculations
		pdfs = pstep(X, params, kvec)
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
	rownames(posterior) = if(isMDF(X)) rownames(X) else names(X)
	
	assignment = apply(posterior, 1, which.max)
	
	## package and return
	
	iteration.stats = Reduce(rbind, iteration.stats)
	rownames(iteration.stats) = NULL
	iteration.stats = as.data.frame(iteration.stats)
	
	retn = namedList(N, D, K, X, npar, npar.hidden, npar.observed, iter, 
		params, stats, weights, pdfs, posterior, assignment, iteration.params, iteration.stats, family, distn, prior, iter.max, dname, dattr, kvec, estep, mstep, pstep, sstep)
	
	class(retn) = "mixmod"
	
	return(retn)
}
