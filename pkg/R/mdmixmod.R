### mdmixmod.R:  define the mdmixmod (multiple data mixture models) function for the lcmix package.

mdmixmod <- function(X, K, K0=min(K), topology=LC_TOPOLOGY, 
	family=NULL, iter.max=LC_ITER_MAX)
{
	## munge data, gather preliminary information, and do sanity checks
		
	if(!is.list(X)) stop("'X' must be a list")
	
	dname = deparse(substitute(X))
	
	if(is.null(names(X)))
		names(X) = paste("X", 1:length(X), sep="")
			# so an unnamed list will now have names "X1", "X2", ...
	dattr = lapply(X, attributes)
	X = lapply(X, function(xelem) drop(as.matrix(xelem)))
		# put all data into matrix or vector form
	
	zvec = names(X)
	Z = length(zvec)
	if(Z < 2) stop("single data source, use mixmod() instead")
	names(zvec) = zvec
		# so the results of lapply() and related functions will always have the names of the elements of X
	
	univariate = !sapply(X, is.matrix)
	ntmp = sapply(zvec, function(z)
		ifelse(univariate[[z]], length(X[[z]]), nrow(X[[z]])))
	if(min(ntmp) != max(ntmp))
		stop("unequal data sizes")
	else
		N = ntmp[[1]]
	
	D = sapply(zvec, function(z) ifelse(univariate[[z]], 1, ncol(X[[z]])))
	
	topology = match.arg(topology)
	
	if(is.null(family))
		family = "normal"
			# we can't use match.arg() here because users may pass in vectors of family names, not just a single name, so we do this instead
	if(length(family) > Z)
		warning("length(family) > Z, truncating to fit")
	if(length(family) != Z)
		family = rep(family, length.out=Z)
	names(family) = zvec
	family = sapply(zvec, function(z)
		match.arg(family[[z]], names(LC_FAMILY)))
		# here is where we finally check for a valid family name
	dummy = lapply(zvec, function(z) .sanity_check_data(X[[z]], family[[z]]))
	
	distn = sapply(zvec, function(z)
		ifelse(univariate[[z]], LC_FAMILY[[(family[[z]])]]$uni, 
			                    LC_FAMILY[[(family[[z]])]]$multi))
	
	if(length(K) > Z)
		warning("length(K) > Z, truncating to fit")
	if(length(K) != Z)
		K = rep(K, length.out=Z)
	names(K) = zvec
	kvec = lapply(zvec, function(z) 1:K[[z]])
	
	k0vec = 1:K0
	
	npar.hidden = switch(topology,
		layered = K0*sum(K-1) + K0-1,
		chained = sum(c(K0, K[1:(Z-1)]) * (K-1)) + K0 - 1
	)
	npar.observed.vec = sapply(zvec, function(z)
		.npar_observed(X[[z]], distn[[z]]))
	npar.observed = sum(K * npar.observed.vec)
	npar = npar.hidden + npar.observed
	bicsub = npar * log(N)
	
	## define internal functions
	
	# weight matrices given PDFs and previous iteration's parameters
	estep <- switch(topology,
	
		layered = function(pdfs, params)
		{
			U = t(t(pdfs$beta) * params$hidden$prob0) / pdfs$gamma
			V = lapply(zvec, function(z) lapply(k0vec, function(k0)
			{
				numerator = t(t(pdfs$G[[z]]) * params$hidden$cprob[[z]][k0,])
				denominator = pdfs$alpha[[z]][,k0]
				U[,k0] * numerator  / denominator
			}))
			W = lapply(zvec, function(z) Reduce("+", V[[z]]))
			
			list(U=U, V=V, W=W)
		},
		
		chained = function(pdfs, params)
		{
			## put parameters and PDFs into mathematical notation
			
			p0 = params$hidden$prob0
			P = params$hidden$probz
			Q = params$hidden$cprob
			R = params$hidden$rprob
			
			alpha = pdfs$alpha
			beta = pdfs$beta
			gamma = pdfs$gamma
			
			## now do the actual calculations
						
			V1 = list(lapply(k0vec, function(k0)
				t(t(beta[[1]])* R[[1]][,k0]) / gamma))
			
			V2plus = lapply(2:Z, function(z) lapply(kvec[[z-1]], function(kzm1)
			{
				term1 = t(alpha[[z-1]][,kzm1] * beta[[z]])
				term2 = Q[[z]][kzm1, ]
				term3 = outer(gamma, P[[z]])
				
				t(term1 * term2) / term3
			}))
			V = c(V1, V2plus)
			names(V) = zvec
			
			U = sapply(V[[1]], rowSums)
			W = lapply(V, function(.) Reduce("+", .))
			
			## send it back
			
			list(U=U, V=V, W=W)
		}
	)
	
	# parameter estimates given weight matrices
	mstep.hidden <- switch(topology,
		
		layered = function(weights)
		{
			usum = colSums(weights$U)
		
			prob0 = usum / N # "p" parameter
			cprob = lapply(zvec, function(z) # "Q" parameter
			{
				numerator = sapply(k0vec, function(k0)
					colSums(weights$V[[z]][[k0]]))
				if(K[[z]] > 1) numerator = t(numerator)
					# K[[z]]=1 produces transposed results for the numerator, hence the above rigamarole instead of just wrapping everything up as "numerator = t(sapply(...))"
				denominator = usum
				as.matrix(numerator / denominator)
					# "as.matrix" to deal with K0=1 case
			})
			
			list(prob0=prob0, cprob=cprob)
		},
		
		chained = function(weights)
		{	
			tH = lapply(weights$V, function(V_z) sapply(V_z, colMeans))
			H = lapply(tH, t)
					
			prob0 = rowSums(H[[1]])
			probz = lapply(H, colSums)
								
			Q1 = list(H[[1]] / prob0)
			Q2plus = lapply(2:Z, function(z) H[[z]] / probz[[z-1]])
			cprob = c(Q1, Q2plus)
			names(cprob) = zvec
			
			rprob = lapply(zvec, function(z) as.matrix(tH[[z]] / probz[[z]]))
				# "as.matrix" to deal with K_z = 1 case; we do not have to do this with cprob because there are calculations in the E-step which require accessing the columns (as such) of the elements of rprob, but no such calculations for the elements of cprob
						
			list(prob0=prob0, probz=probz, cprob=cprob, rprob=rprob)
		}
	)
	mstep.observed <- lapply(zvec, function(z) .mstep_observed(distn[[z]]))
	mstep <- function(weights)
	{
		hidden = mstep.hidden(weights)
		observed = lapply(zvec, function(z)
			mstep.observed[[z]](X[[z]], kvec[[z]], weights$W[[z]]))
		
		list(hidden=hidden, observed=observed)
	}
	
	# PDFs given parameter estimates
	pstep.main <- switch(topology,
		
		layered = function(params, G)
		{
			alpha = lapply(zvec, function(z)
				G[[z]] %*% t(params$hidden$cprob[[z]]))
			beta = Reduce("*", alpha)
			gamma = drop(beta %*% params$hidden$prob0) # = f(X)
				# "drop" so we have an N-length vector rather than a N-by-1 matrix
			
			list(alpha=alpha, beta=beta, gamma=gamma)
		},
		
		chained = function(params, G)
		{
			## put parameters into mathematical notation
			
			p0 = params$hidden$prob0
			P = params$hidden$probz
			Q = params$hidden$cprob
			R = params$hidden$rprob
			
			## now do the actual calculations
			
			alpha = vector(mode="list", length=Z)
			alpha[[1]] = t(t(G[[1]]) * P[[1]])
			for(z in 2:Z)
				alpha[[z]] = G[[z]] * (alpha[[z-1]] %*% Q[[z]])
			names(alpha) = zvec
			
			beta = vector(mode="list", length=Z)
			beta[[Z]] = t(t(G[[Z]]) * P[[Z]])
			for(z in (Z-1):1)
				beta[[z]] = G[[z]] * (beta[[z+1]] %*% R[[z+1]])
			names(beta) = zvec
			
			gamma = rowSums(alpha[[Z]])
						
			## send it back
			
			list(alpha=alpha, beta=beta, gamma=gamma)
		}
	)
	pstep.observed <- lapply(zvec, function(z) .pstep_observed(distn[[z]]))
	pstep <- function(params)
	{
		
		G = lapply(zvec, function(z)
		{
			res = pstep.observed[[z]](X[[z]], kvec[[z]], params$observed[[z]])
			pmax(res, LC_EPSEPS) # see mstep() in mixmod() to explain this
		})
		
		c(list(G=G), pstep.main(params, G))
	}
	
	# iteration statistics given iteration and results of previous steps
	sstep <- function(iter, weights, params, pdfs)
	{
		llik = sum(log(pdfs$gamma))
		qval = sum(.qfun_mdmixmod(weights, params, pdfs))
		bic = 2*llik - bicsub
		
		c(iter=iter, llik=llik, qval=qval, bic=bic)
	}
	
	## initialize
	
	# quasi E-step
	Xinit = lapply(zvec, function(z)
		if(family[[z]] %in% LC_NONNEGFAM) log(X[[z]]) else X[[z]])
			# convert to (-Inf, Inf) for Gaussian weights
	W = lapply(zvec, function(z) .init_weights(Xinit[[z]], K[[z]]))
	if(topology == "layered") { # give equal weight to all data sources
		Uinit = lapply(Xinit, function(Xinit_z) .init_weights(Xinit_z, K0))
		U = Reduce("*", Uinit)
		U = U / rowSums(U)
	} else if(topology == "chained") { # calculate from first data source
		U = .init_weights(Xinit[[1]], K0)
	} # "else":  bad topology, to be caught by match.arg() above
	V = switch(topology,
			layered = lapply(zvec, function(z)
						lapply(k0vec, function(k0) U[,k0] * W[[z]])),
			chained = c(list(lapply(k0vec, function(k0) U[,k0] * W[[1]])),
			            lapply(2:Z, function(z)
			            	lapply(kvec[[z-1]], function(kzm1)
			            		W[[z-1]][,kzm1] * W[[z]])))
			           
	)
	names(V) = zvec
		# the "chained" choice in the switch statement above requires this, and it will have no effect on the "layered" choice, so it's easier to just put it outside the switch statement
	weights = list(U=U, V=V, W=W)
	
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
	
	posterior = weights$U
	
	assignment = apply(posterior, 1, which.max)
	
	## package and return
	
	iteration.stats = Reduce(rbind, iteration.stats)
	rownames(iteration.stats) = NULL
	iteration.stats = as.data.frame(iteration.stats)
	
	retn = list(N                = N,
	            Z                = Z,
	            D                = D,
	            K                = K,
	            K0               = K0,
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
	            topology         = topology,
	            family           = family,
	            distn            = distn,
	            iter.max         = iter.max,
	            dname            = dname,
	            dattr            = dattr)
	class(retn) = c("mdmixmod", "mixmod")
	
	return(retn)
}