### invisible_utilities:  invisible utility functions for the lcmix package

# Replace bad numeric values (NA, NaN, Inf) with 0
.bad2zero <- function(x)
{
	x[is.na(x)] = 0
	x[is.nan(x)] = 0
	x[is.infinite(x)] = 0
	
	return(x)
}

# Model E-step functions:  given no arguments; return functions of parameters and PDF matrices in the single-data case; or functions of parameters, PDF matrices, vector of data source names, and vectors of possible component values in the multiple-data case; which calculate weights.
.estep_mixmod <- function()
{
	function(pdfs, params)
	{
		lW = as.matrix(t(t(pdfs$lG) + params$hidden$lprob) - pdfs$lgamma)
			# "as.matrix" to deal with the pathological K=1 case
		W = exp(lW)
		
		namedList(W)
	}
}
.estep_mdmixmod_layered <- function()
{
	function(pdfs, params, zvec, kvec, k0vec)
	{	
		## put parameters and PDFs into mathematical notation
		
		lp0 = params$hidden$lprob0
		lQ = params$hidden$lcprob
		
		lG = pdfs$lG
		lalpha = pdfs$lalpha
		lbeta = pdfs$lbeta
		lgamma = pdfs$lgamma
		
		## now do the actual calculations
		
		lU = t(t(lbeta) + lp0) - lgamma
		U = exp(lU)
		
		lV = lapply(zvec, function(z) lapply(k0vec, function(k0)
		{
			lnumerator = t(t(lG[[z]]) + lQ[[z]][k0,])
			ldenominator = lalpha[[z]][,k0]
			lU[,k0] + lnumerator - ldenominator
		}))
		V = rlapply(lV, exp, "matrix")
		
		W = lapply(V, listSum)
		
		## send it back

		namedList(U, V, W)
	}
}
.estep_mdmixmod_chained <- function()
{
	function(pdfs, params, zvec, kvec, k0vec)
	{
		## put parameters and PDFs into mathematical notation
		
		lp0 = params$hidden$lprob0
		lP = params$hidden$lprobz
		lQ = params$hidden$lcprob
		lR = params$hidden$lrprob
		
		lalpha = pdfs$lalpha
		lbeta = pdfs$lbeta
		lgamma = pdfs$lgamma
		
		## now do the actual calculations
		
		lfy1gx = lbeta[[1]] - lgamma
		tlfy1gx = t(lfy1gx)
		lV1 = lapply(k0vec, function(k0) t(tlfy1gx + lR[[1]][,k0]))
		
		lV2plus = lapply(2:length(zvec), function(z)
		{
			term3 = outer(lgamma, lP[[z]], "+") # more efficient to do this here
			
			lapply(kvec[[z-1]], function(kzm1)
			{ 
				term1 = t(lalpha[[z-1]][,kzm1] + lbeta[[z]])
				term2 = lQ[[z]][kzm1, ]
				
				t(term1 + term2) - term3
			})
		})
		
		#lV426:  stay away!
			
		lV = c(list(lV1), lV2plus)
		names(lV) = zvec
		V = rlapply(lV, exp, "matrix")
		
		U = sapply(V[[1]], rowSums)
		W = lapply(V, listSum)
		
		## send it back
		
		namedList(U, V, W)
	}
}

# Return a fully formed "dvec" vector given a matrix or data frame
.extract_dvec <- function(X)
{ 
	retn = 1:ncol(X)
	names(retn) = colnames(X)
	
	return(retn)
}

# Quickly calculates the inverse of double (NOT INTEGER) valued matrix X; does no sanity checking, offers no options, and in general should only be used when you're sure this is what you want.  To iterate over a list, and really speed things up, you can pass in the identity matrix; e.g, if Y is a list of 3x3 matrices, you can write "I3 = diag(3) ; res = lapply(Y, .fast_inv, I3)" which will be considerably faster than just "res = lapply(Y, .fast_inv)".  The fastest way to convert X from integer to double is just to call ".fast_inv(X+0.0)".  The ".list_fast_inv" wrapper is for convenience and speed; applied to a list, it returns a list.  Finally, "*sym()" guarantees that the return value will be symmetric, because inversion can sometimes introduce slight asymmetries; these methods should of course be used only with matrices that are supposed to be symmetric, like covariance matrices.
#*# TO DO:  either find a way to write equally fast code using pure R (no calls to "La_dgesv" or other base package ".Call" functions) or link directly to the relevant LAPACK functions.  Per <http://stats.stackexchange.com/questions/14951/efficient-calculation-of-matrix-inverse-in-r>, "chol2inv(chol(.))" may be much faster than "solve(.)" (and guarantee symmetry of the results?)
#*# .fast_inv <- function(X, ID=diag(nrow(X)))
#*# 	.Call("La_dgesv", X, ID, LC_EPS, PACKAGE="base")
.fast_inv  <- function(X, ID=diag(nrow(X))) solve(X, ID)
.fast_inv_sym <- function(X, ID=diag(nrow(X)))
{
	res = .fast_inv(X, ID)
	(res + t(res)) / 2.0
}
.list_fast_inv <- function(L)
	lapply(L, .fast_inv, ID=diag(nrow(L[[1]])))
.list_fast_inv_sym <- function(L)
	lapply(L, .fast_inv_sym, ID=diag(nrow(L[[1]])))

# Quickly calculates the logarithm of the absolute value of the determinant of double (NOT INTEGER) valued matrix X; does no sanity checking, offers no options, and in general should only be used when you're sure this is what you want.  The ".vec_fast_labs_det" wrapper is for convenience; applied to a list, it returns a vector.
#*# TO DO:  either find a way to write equally fast code using pure R (no calls to "det_ge_real" or other base package ".Call" functions) or link directly to the relevant LAPACK functions.
#*# .fast_labs_det <- function(X)
#*# 	.Call("det_ge_real", X, TRUE, PACKAGE="base")$modulus
.fast_labs_det <- function(X) determinant(X, logarithm=TRUE)$modulus
.vec_fast_labs_det <- function(L) sapply(L, .fast_labs_det)

# Calculate the initial weight matrix ("W") for data "X" with "K" means, means being chosen appropriately. The matrix that is returned is N-by-K.
.init_weights <- function(X, K, prior=NULL, rowsmask=NULL)
{
	X = apply(as.matrix(X), 2, rank)
	mu = uniformQuantiles(X, K, decreasing=TRUE)
	
	D = sapply(1:K, function(k) manhattan(X, mu[k,])) # distance matrix
	S = 1/(D+LC_EPS) # similarity matrix
	if(!is.null(prior)) S = t(t(S) * prior)
	
	retn = S / rowSums(S)
	if(!is.null(rowsmask)) retn[rowsmask$rows,] = rowsmask$mask
	return(retn)
}

# TESTING

.inititer_norm <- function(x, mu)
{
	# initialization (partial M- and P-steps)
	sigma = sqrt(colMeans(outer(x, mu, "-")^2))
		# we do E[X - mu]^2 here to avoid negative variances since the "means" are randomly chosen; below in the M-step, we use the simpler E[X^2] - mu^2 method
	G = mapply(function(mk, sk) dnorm(x, mk, sk), mu, sigma)
	
	# E- and M-steps
	W = rowDiv(G)
	csw = colSums(W)
	p = colMeans(W)
	mu = drop(x %*% W) / csw
	sigma = sqrt(colSums(x^2 * W)/csw - mu^2)
	
	# P- and S-steps
	G = mapply(function(mk, sk) dnorm(x, mk, sk), mu, sigma)
	gamma = drop(G %*% p)
	lgamma = log(gamma)
	llik = sum(lgamma)
	
	# prevent inversion
	muorder = order(mu, decreasing=TRUE)
	W = W[,muorder]
	mu = mu[muorder]
	p = p[muorder]
	
	# entropy
	entr = -2 * sum(W*log(W))
	
	# package and return
	namedList(W, mu, p, llik, entr)
}

.inititer_mvnorm <- function(X, mu)
{
	# gather useful statistics and munge data
	N = nrow(X)
	K = nrow(mu)
	kvec = 1:K
	tX = t(X)
	
	# initialization (partial M- and P-steps)
	Sigma = lapply(kvec, function(k) tcrossprod(tX - mu[k,]) / N)
		# see note on the initial "sigma" in .inititer_norm()
	G = sapply(kvec, function(k) dmvnorm(X, mu[k,], Sigma[[k]]))
	
	# E- and M-steps
	W = rowDiv(G)
	csw = colSums(W)
	p = colMeans(W)
	mu = t(tX %*% W) / csw
	Sigma = lapply(kvec, function(k)
		crossprod(X, W[,k]*X)/csw[k] - tcrossprod(mu[k,]))
	
	# P- and S-steps
	G = sapply(kvec, function(k) dmvnorm(X, mu[k,], Sigma[[k]]))
	gamma = drop(G %*% p)
	lgamma = log(gamma)
	llik = sum(lgamma)
	
	# prevent inversion
	muorder = order(rowSums(mu), decreasing=TRUE)
	W = W[,muorder]
	mu = mu[muorder,]
	p = p[muorder]
	
	# entropy
	entr = -2 * sum(W*log(W))
	
	# package and return
	namedList(W, mu, p, llik, entr)
}

.init_weights2 <- function(X, K,  nstart=100, seed=123) # nstart=100, seed=NULL
{
	X = drop(as.matrix(X))
	
	multivariate = is.matrix(X)
	iterfun <- if(multivariate) .inititer_mvnorm else .inititer_norm
	sampfun <- if(multivariate) function()
		rowSort(rowSample(X, K), decreasing=TRUE)
	else function()
		sort(sample(X, K), decreasing=TRUE)
	
	N = ifelse(multivariate, nrow(X), length(X))
	if(multivariate) tX = t(X)
	if(!is.null(seed)) set.seed(seed)
	
	res = lapply(1:nstart, function(.) iterfun(X, sampfun()))
	llik = sapply(res, function(r) r$llik)
	idx = which.max(llik)
	
	retn = res[[idx]]$W
	attr(retn, "mu") = res[[idx]]$mu
	attr(retn, "p") = res[[idx]]$p
	attr(retn, "llik") = res[[idx]]$llik
	return(retn)
}

# /TESTING

# Given weight matrices "A" and "B", with column means a and b, calculate the transition matrix C such that:
# (1) a C = b
# (2) the value sum((A C - B)^2) is minimized
# using .dependentDiscreteCpdf (implying positive dependence).
.least_squares_transition_matrix <- function(A, B)
{
	a = colMeans(A)
	b = colMeans(B)
	toptim <- function(d) sum(((A %*% .dependentDiscreteCpdf(a, b, d)) - B)^2)
	
	.dependentDiscreteCpdf(a, b, optimize(toptim, c(0,1))$minimum)
}

# Deal with outliers in log-pdf-by-component matrices so extreme their density is effectively 0 for all components; the value -300 is sufficient to prevent instability but not actually change the results in any meaningful way
.lGfilter <- function(lG)
{
	if(min(lG) < -300) {
		K = ncol(lG)
		lG[(rowSums(lG < -300) == K), ] = -300
	}
	return(lG)
}
	
# Workhorse function for (unweighted) MLE of covariance; takes a D-by-N matrix, a D-length vector of means, and N, and returns a MLE covariance matrix.
.mleTcov <- function(tX, mu, N) { tXadj = tX - mu ; tcrossprod(tXadj) / N }

# Model M-step functions:  given distribution(s), return functions of data and weights (and data source names in the multiple-data-model case) which calculate parameters.
.mstep_mixmod <- function(distn)
{
	mstep.observed <- .mstep_observed(distn)
	
	function(X, weights, theta=NULL, prior=NULL) {
		
		prob = if(is.null(prior)) colMeans(weights$W) else(prior)
		lprob = log(prob)
		hidden = namedList(prob, lprob)
		
		observed = mstep.observed(X, weights$W, theta=theta$observed)
		
		namedList(hidden, observed)
	}
}
.mstep_mdmixmod_layered <- function(distn)
{
	mstep.observed <- mapply(.mstep_observed, distn)
	
	function(X, weights, K, kvec, K0, k0vec, zvec, params=NULL, prior=NULL)
	{
		prob. = colMeans(weights$U)
		
		prob0 = if(is.null(prior)) prob. else(prior)
		lprob0 = log(prob0)
		lcprob = lapply(weights$V, function(Vz)
			log(t(sapply(Vz, colMeans))) - log(prob.))
		cprob = lapply(lcprob, exp)
		hidden = namedList(prob0, cprob, lprob0, lcprob)
		
		observed = lapply(zvec, function(z)
			mstep.observed[[z]](X[[z]], weights$W[[z]],
				theta=params$observed[[z]]))
		namedList(hidden, observed)
	}
}
.mstep_mdmixmod_chained <- function(distn)
{
	mstep.observed <- mapply(.mstep_observed, distn)
	
	function(X, weights, K, kvec, K0, k0vec, zvec, params=NULL, prior=NULL)
	{
		Z = length(zvec)
	
		tH = lapply(weights$V, function(V_z) sapply(V_z, colMeans))
		ltH = lapply(tH, log)
		H = lapply(tH, t)
		lH = lapply(ltH, t)
		
		prob. = rowSums(H[[1]])
		prob0 = if(is.null(prior)) prob. else prior
		lprob0 = log(prob0)
		probz = lapply(H, colSums)
		lprobz = lapply(probz, log)
							
		lQ1 = lH[[1]] - log(prob.)
		lQ2plus = lapply(2:Z, function(z) lH[[z]] - lprobz[[z-1]])
		lcprob = c(list(lQ1), lQ2plus)
		names(lcprob) = zvec
		cprob = lapply(lcprob, exp)
		
		lrprob = lapply(zvec, function(z) ltH[[z]] - lprobz[[z]])
		rprob = lapply(lrprob, exp)
					
		hidden = namedList(prob0, probz, cprob, rprob,
		     lprob0, lprobz, lcprob, lrprob)
				
		observed = lapply(zvec, function(z)
			mstep.observed[[z]](X[[z]], weights$W[[z]],
				theta=params$observed[[z]]))
		
		namedList(hidden, observed)
	}
}

# Return a function which performs weighted observed-distribution parameter estimation for the single-data mixture model M-step (can also be used in the multiple-data mixture model if handled appropriately.)
.mstep_observed <- function(distn) get(sprintf(".thetahat_%s", distn))

# Combined E- and M-steps for univariate Pearson Type VII parameter estimation.  Takes data, transposed data, weights, and a list "theta" containing the sum of weights, log of sum of weights, and current parameter values.  Returns updated theta with an additional element:  the value of the objective function.
.mvpvii_naive_emstep <- function(X, tX, w, sw, theta)
{
	# E-step
	
	D = ncol(X)
	N = nrow(X)
	X.prime = t(t(X) - theta$mean)
	iSigma = .fast_inv_sym(theta$scale)
	I = diag(D)
	
	alpha.prime = theta$shape + 0.5
	Psi = lapply(1:N, function(n) tcrossprod(X.prime[n,]))
	Lambda.prime = lapply(Psi, function(Psin) I + 0.5*Psin%*%iSigma)
	
	iLambda.prime = .list_fast_inv_sym(Lambda.prime)
	T.prime = lapply(iLambda.prime, "*", y=alpha.prime)
	mgammaterm = sum(digamma(alpha.prime - (0:(D-1))/2))
	labsdetterm = .vec_fast_labs_det(Psi)
	t.star = mgammaterm - labsdetterm
	
	wtp = lapply(1:N, function(n) w[n] * T.prime[[n]])
	swtp = listSum(wtp)
	wtpx = lapply(1:N, function(n) wtp[[n]] %*% X[n,])
	swtpx = listSum(wtpx)
	wpsitp = lapply(1:N, function(n) w[n] * Psi[[n]] %*% T.prime[[n]])
	swpsitp = listSum(wpsitp)
	swts = sum(w * t.star)

print(swpsitp)
	
	# M-step
		
	mu = drop(.fast_inv_sym(swtp) %*% swtpx)
	Sigma = swpsitp / sw
	
	toptim <- function(alpha) -lmgamma(alpha, D)*sw - alpha*swts
	res = nlminb(theta$shape, toptim, lower=LC_EPS+(D-1)/2)
	alpha = res$par
	
	# send it back
	list(mean=mu, scale=Sigma, shape=alpha, obj=res$objective)
}
.mvpvii_emstep <- function(X, tX, w, sw, theta)
{
	# E-step
	
	shape.prime = theta$shape + 0.5*ncol(X)
	tX.prime = tX - theta$mean
	X.prime = t(tX.prime)
	delta = rowSums((X.prime %*% .fast_inv_sym(theta$scale)) * X.prime)
	rate.prime = 1 + 0.5*delta
	
	tprime  = shape.prime / rate.prime
	tstar   = digamma(shape.prime) - log(rate.prime)
	
	# M-step
	
	wtp = w * tprime
	mu = colSums(wtp * X) / sum(wtp)
	
	Sigma = tX.prime %*% (wtp * X.prime) / sw
	
	swts = sum(w * tstar)
	toptim <- function(alpha) lgamma(alpha)*sw - alpha*swts
	res = nlminb(theta$shape, toptim, lower=1+LC_EPS)
	alpha = res$par
	
	# send it back
	list(mean=mu, scale=Sigma, shape=alpha, obj=res$objective)
}

# Initial parameter estimation for multivariate Pearson Type VII.  Takes data, weights.  Returns the estimated parameters.
.mvpvii_init_params <- function(X, w)
{
	theta.mvnorm = thetahat.mvnorm(X, w)

	mean = theta.mvnorm$mean
	shape = 2
	scale = theta.mvnorm$cov * shape
	
	namedList(mean, scale, shape)
}

# Calculate the number of parameters for distribution of observed data in a single-data mixture model
.npar_observed <- function(distn, K, D)
{
	K * switch(distn,
		norm    = 2,
		mvnorm  = D + D*(D+1)/2,
		weisd   = 2,
		mvweisd = 2*D + D*(D-1)/2,
		gamma   = 2,
		mvgamma = 2*D + D*(D-1)/2,
		exp     = 1,
		mvexp   = D + D*(D-1)/2,
		pvii    = 3,
		mvpvii  = D + D*(D+1)/2 + 1,
		stop(sprintf("unknown distribution '%s'", distn))
	)
}

# Return a nicely formatted string version of a family name, suitable for printing
.pretty_family <- function(family)
{
	switch(family,
		weibull = "Weibull",
		pvii = "PVII",
		family)
}

# Model P-step functions:  given distribution(s), return functions of data, number(s) of components, and parameters (and data source names in the multiple-data-model case) which calculate PDF matrices
.pstep_mixmod <- function(distn)
{
	pstep.observed <- .pstep_observed(distn)
	
	function(X, params, kvec)
	{
		lG = .lGfilter(pstep.observed(X, kvec, params$observed))
		G = exp(lG) # g_n,k = f(x_n | y=k)
		gamma = drop(G %*% params$hidden$prob) # f(X)
			# "drop" so we have an N-length vector rather than a N-by-1 matrix
		lgamma = log(gamma)
		
		namedList(G, gamma, lG, lgamma)
	}
}
.pstep_mdmixmod_layered <- function(distn)
{
	pstep.observed <- mapply(.pstep_observed, distn)
	
	function(X, params, kvec, zvec)
	{
		lG = lapply(zvec, function(z) .lGfilter(pstep.observed[[z]](
			X[[z]], kvec[[z]], params$observed[[z]])))
		G = lapply(lG, exp)
		
		alpha = lapply(zvec, function(z)
			G[[z]] %*% t(params$hidden$cprob[[z]]))
		lalpha = lapply(alpha, log)
		lbeta = listSum(lalpha)
		beta = exp(lbeta)
		gamma = drop(beta %*% params$hidden$prob0) # = f(X)
			# "drop" so we have an N-length vector rather than a N-by-1 matrix
		lgamma = log(gamma)
		
		namedList(G, lG, alpha, beta, gamma, lalpha, lbeta, lgamma)
	}
}
.pstep_mdmixmod_chained <- function(distn)
{
	
	pstep.observed <- mapply(.pstep_observed, distn)
	
	function(X, params, kvec, zvec)
	{	
		## initial component-specific PDF calculations
		
		lG = lapply(zvec, function(z) .lGfilter(pstep.observed[[z]](
			X[[z]], kvec[[z]], params$observed[[z]])))
		
		G = lapply(lG, exp)
		
		## put parameters into mathematical notation
			
		P = params$hidden$probz
		Q = params$hidden$cprob
		R = params$hidden$rprob
		
		lP = params$hidden$lprobz
		
		## now do the actual calculations
		
		Z = length(zvec)
		
		alpha = vector(mode="list", length=Z)
		lalpha = vector(mode="list", length=Z)
		lalpha[[1]] = t(t(lG[[1]]) + lP[[1]])
		alpha[[1]] = exp(lalpha[[1]])
		for(z in 2:Z)
		{
			lalpha[[z]] = lG[[z]] + log(alpha[[z-1]] %*% Q[[z]])
			alpha[[z]] = exp(lalpha[[z]])
		}
		names(alpha) = zvec
		
		beta = vector(mode="list", length=Z)
		lbeta = vector(mode="list", length=Z)
		lbeta[[Z]] = t(t(lG[[Z]]) + lP[[Z]])
		beta[[Z]] = exp(lbeta[[Z]])
		for(z in (Z-1):1)
		{
			lbeta[[z]] = lG[[z]] + log(beta[[z+1]] %*% R[[z+1]])
			beta[[z]] = exp(lbeta[[z]])
		}
		names(beta) = zvec
		
		gamma = rowSums(alpha[[Z]])
		lgamma = log(gamma)
		
		## send it back
		
		namedList(G, lG, alpha, beta, gamma, lalpha, lbeta, lgamma)
	}
}

# Return a function which performs PDF estimation for the observed data given the components of a single-data mixture model (can also be used in multiple-data mixture models if handled appropriately.)
.pstep_observed <- function(distn)
{
	switch(distn,
		
		norm = function(X, kvec, obspar) 
			sapply(kvec, function(k)
				dnorm(X, obspar$mean[k], obspar$sd[k], log=TRUE)),
		mvnorm = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvnorm(X, obspar$mean[k,], obspar$cov[[k]], log=TRUE)),
		
		weisd = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dweisd(X, obspar$shape[k], obspar$decay[k], log=TRUE)),
		mvweisd = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvweisd(X, obspar$shape[k,], obspar$decay[k,],
					obspar$corr[[k]], log=TRUE)),
		
		gamma = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dgamma(X, obspar$shape[k], obspar$rate[k], log=TRUE)),
		mvgamma = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvgamma(X, obspar$shape[k,], obspar$rate[k,],
					obspar$corr[[k]], log=TRUE)),
		
		exp = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dexp(X, obspar$rate[k], log=TRUE)),
		mvexp = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvexp(X, obspar$rate[k,],
					obspar$corr[[k]], log=TRUE)),
		
		pvii = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dpvii(X, obspar$mean[k], obspar$scale[k],
					obspar$shape[k], log=TRUE)),
		mvpvii = function(X, kvec, obspar)
			sapply(kvec, function(k)
				dmvpvii(X, obspar$mean[k,], obspar$scale[[k]],
					obspar$shape[k], log=TRUE))
	
	)
}

# Combined E- and M-steps for univariate Pearson Type VII parameter estimation.  Takes data, weights, and a list "theta" containing the sum of weights, log of sum of weights, and current parameter values.  Returns updated theta with an additional element:  the value of the objective function.
.pvii_emstep <- function(x, w, sw, theta)
{	
	# E-step
	
	shape.prime = theta$shape + 0.5
	Delta = (x - theta$mean)^2
	delta = Delta / theta$scale
	rate.prime  = 1 + 0.5*delta
	
	tprime  = shape.prime / rate.prime
	tstar   = digamma(shape.prime) - log(rate.prime)
	
	# M-step
	
	wtp = w * tprime
	mu = sum(wtp * x) / sum(wtp)
	
	sigmasq = sum(wtp * Delta) / sw
	
	swts = sum(w * tstar)
	toptim <- function(alpha) lgamma(alpha)*sw - alpha*swts
	res = nlminb(theta$shape, toptim, lower=1+LC_EPS)
	alpha = res$par
	
	# send it back
	
	list(mean=mu, scale=sigmasq, shape=alpha, obj=res$objective)
}

# Initial parameter estimation for univariate Pearson Type VII.  Takes data, weights.  Returns the estimated parameters.
.pvii_init_params <- function(x, w)
{
	theta.norm = thetahat.norm(x, w)
	
	mean = theta.norm$mean
	shape = 2
	scale = theta.norm$var * shape
	
	namedList(mean, scale, shape)
}

# "Q-function" (expected log-likelihood) for single-data (.qfun_mixmod()) and multiple-data (.qfun_mdmixmod()) mixture models. If map==FALSE, the expected value will be returned; if TRUE, the MAP (maximum a posteriori) estimate will be returned.  The .bad2zero() calls are to deal with what happens when we try to take the log of 0 probabilities.
.qfun_mixmod <- function(weights, params, pdfs, map=FALSE)
{
	W = weights$W
	if(map) W = isRowMax(W)
	
	hidden = sum(.bad2zero(W %*% params$hidden$lprob))
	observed = sum(.bad2zero(weights$W * pdfs$lG))
	
	namedVector(hidden, observed)
}
.qfun_mdmixmod <- function(weights, params, pdfs, topology, map=FALSE)
{
	U = weights$U
	V = weights$V
	W = weights$W
	
	zvec = names(W)
	names(zvec) = zvec
	Z = length(zvec)
	
	if(map) {
		U = isRowMax(U)
		V = lapply(V, .v2map)
		W = lapply(W, isRowMax)
	}
	
	lp = params$hidden$lprob0
	lQ = params$hidden$lcprob
	lG = pdfs$lG
	
	term1 = sum(.bad2zero(U %*% lp))
	
	term2 = sum(sapply(zvec, function(z)
		sum(sapply(1:length(V[[z]]), function(kparens)
			sum(.bad2zero(V[[z]][[kparens]] %*% lQ[[z]][kparens,]))))))
	
	term3 = sum(sapply(zvec, function(z)
		sum(.bad2zero(W[[z]] * lG[[z]]))))
	
	c(hidden=term1+term2, observed=term3)
}

# Weighted multivariate standard normal (means all equal to 0, variances all equal to 1) correlation (equivalent to covariance) matrix estimation.
.rhohat <- function(X, w) cov2cor(t(X) %*% (w*X) / sum(w))
	# the columns of X are already assumed to have mean 0, so no need to estimate means and center; "cov2cor" is a hack to make sure the result is a valid correlation matrix, and fixing it is TO DO

# Return a list containing the outer products of the rows of matrix X.
.tcrossprodcols <- function(X) lapply(1:ncol(X), function(n) tcrossprod(X[,n]))
.tcrossprodrows <- function(X) lapply(1:nrow(X), function(n) tcrossprod(X[n,]))

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

# Turn a weight matrix into a matrix of almost-0 and almost-1 values, where the maximum in each row gets assigned almost-1 and all other positions almost-0
.sharpen_weight_matrix <- function(W)
{
	res = diag(ncol(W))[apply(W, 1, which.max),] + LC_EPS
	
	res / rowSums(res)
}

# Hidden data simulation:  take a number of samples to simulate, a distribution, and hidden parameters (and a topology for .simulate_hidden_mdmixdata()) and return hidden data (Y for .simulate_hidden_mixdata(), a list with elements $Y and $Y0 for .simulate_hidden_mdmixdata()).  No sanity checking is done, so everything has to be properly prepared; e.g., see calls to .simulate_hidden_mdmixdata() from simulateMdmixdata() to see how arguments should be set up.
.simulate_hidden_mixdata <- function(n, par)
{
	rcategorical(n, par$prob)
}
.simulate_hidden_mdmixdata <- function(n, par, topology)
{
	# gather information
	K0 = length(par$prob0)
	K = sapply(par$cprob, ncol)
	Z = length(par$cprob)
	zvec = names(par$cprob)
	names(zvec) = zvec

	# simulate top-level
	Y0 = rcategorical(n, par$prob0)
	
	# simulate lower-level
	if(topology == "layered") {
	
		Y = lapply(zvec, function(z)
		{
			res = vector(mode="numeric", length=n)
			for(k0 in 1:K0)
			{
				isk0 = (Y0 == k0)
				nisk0 = sum(isk0)
				res[isk0] = rcategorical(nisk0, par$cprob[[z]][k0,])
			}
			return(res)
		})
	
	} else if(topology == "chained") {
	
		Y = vector(mode="list", length=Z)
		Y[[1]] = vector(mode="numeric", length=n)
		
		for(k0 in 1:K0)
		{
			isk0 = (Y0 == k0)
			nisk0 = sum(isk0)
			Y[[1]][isk0] = rcategorical(nisk0, 
				par$cprob[[1]][k0,])
		}
		
		if(Z > 1) for(z in 2:Z)
		{
			Y[[z]] = vector(mode="numeric", length=n)
			for(kzm1 in 1:K[[z-1]])
			{
				iskzm1 = (Y[[z-1]] == kzm1)
				niskzm1 = sum(iskzm1)
				Y[[z]][iskzm1] = rcategorical(niskzm1, 
					par$cprob[[z]][kzm1,])
			}
		}
		
		names(Y) = zvec
	
	} else stop(sprintf("unrecognized topology '%s'", topology))
	
	# package it up, send it back
	namedList(Y, Y0)
}

# Model S-step functions:  given no arguments, return functions of iteration, weights, parameters, PDFs, and BIC adjustment factor (size of parameter space times log of sample size) which calculate iteration statistics
.sstep_mixmod <- function(bicsub)
{
	function(iter, weights, params, pdfs)
	{
		llik = sum(pdfs$lgamma)
		qval = sum(.qfun_mixmod(weights, params, pdfs))
		bic = 2*llik - bicsub
		iclbic = 2*qval - bicsub
		
		namedVector(iter, llik, qval, bic, iclbic)
	}
}
.sstep_mdmixmod <- function()
{
	function(iter, weights, params, pdfs, bicsub)
	{
		llik = sum(pdfs$lgamma)
		qval = sum(.qfun_mdmixmod(weights, params, pdfs))
		bic = 2*llik - bicsub
		iclbic = 2*qval - bicsub
		
		namedVector(iter, llik, qval, bic, iclbic)
	}
}

# Weighted parameter estimation functions for mixtures.  If argument "theta" is not NULL, it should be a list of parameters with the same names as those to be estimated, reflecting previous iteration's values of those parameters.  (Note that it will have no effect in the case of distributions with closed form MLEs.)
.thetahat_exp <- function(x, W, aslist=TRUE, theta=NULL)
{
	res = sapply(1:ncol(W), function(k)
		thetahat.exp(x, W[,k], aslist=FALSE))
	
	lambda = as.vector(res)
		# the kth element of lambda is the rate for the kth component
	
	if(aslist)
		list(rate=lambda)
	else
		c(rate=lambda)
}
.thetahat_mvexp <- function(X, W, aslist=TRUE, theta=NULL)
{
	kvec = 1:ncol(W)
	
	res = sapply(kvec, function(k)
		thetahat.mvexp(X, W[,k], aslist=TRUE), simplify=FALSE)
	
	lambda = t(sapply(kvec, function(k) res[[k]]$rate))
		# the kth row of lambda is the vector of the column rates of X for the kth component
	rho = lapply(kvec, function(k) res[[k]]$corr)
	
	if(aslist)
		list(rate=lambda, corr=rho)
	else
		c(rate=lambda, corr=rho)
}
.thetahat_gamma <- function(x, W, aslist=TRUE, theta=NULL)
{
	res = sapply(1:ncol(W), function(k)
		thetahat.gamma(x, W[,k], aslist=FALSE, shape.guess=theta$shape[k]))
	
	sigma = res["shape",]
	lambda = res["rate",]
		# the kth element of sigma is the shape for the kth component, and similarly for lambda (rate)
	
	if(aslist)
		list(shape=sigma, rate=lambda)
	else
		c(shape=sigma, rate=lambda)
}
.thetahat_mvgamma <- function(X, W, aslist=TRUE, theta=NULL)
{
	kvec = 1:ncol(W)
	
	res = sapply(kvec,
	             function(k)
	             	thetahat.mvgamma(X, W[,k], aslist=TRUE, 
	             		shape.guess=theta$shape[k,]),
		         simplify=FALSE)
	
	sigma = t(sapply(kvec, function(k) res[[k]]$shape))
		# the kth row of sigma is the vector of the column shapes of X for the kth component
	lambda = t(sapply(kvec, function(k) res[[k]]$rate))
		# the kth row of lambda is the vector of the column rates of X for the kth component
	rho = lapply(kvec, function(k) res[[k]]$corr)
	
	if(aslist)
		list(shape=sigma, rate=lambda, corr=rho)
	else
		c(shape=sigma, rate=lambda, corr=rho)
}
.thetahat_norm <- function(x, W, aslist=TRUE, theta=NULL)
{
	res = sapply(1:ncol(W), function(k)
		thetahat.norm(x, W[,k], aslist=FALSE))
	
	mu = res["mean",]
	sigmasq = res["var",]
	sigma = res["sd",]
		# the kth element of mu is the mean for the kth component, and similarly for sigmasq (variance) and sigma (std. dev.)
	
	if(aslist)
		list(mean=mu, var=sigmasq, sd=sigma)
	else
		c(mean=mu, var=sigmasq, sd=sigma)
}
.thetahat_mvnorm <- function(X, W, aslist=TRUE, theta=NULL)
{
	
	kvec = 1:ncol(W)
	
	res = sapply(kvec, function(k)
		thetahat.mvnorm(X, W[,k], aslist=TRUE), simplify=FALSE)
	
	mu = t(sapply(kvec, function(k) res[[k]]$mean))
		# the kth row of mu is the vector of the column means of X for the kth component
	Sigma = lapply(kvec, function(k) res[[k]]$cov)
		# the kth list element of Sigma is the covariance matrix for the kth component
		
	if(aslist)
		list(mean=mu, cov=Sigma)
	else
		c(mean=mu, cov=Sigma)
}
.thetahat_pvii <- function(x, W, aslist=TRUE, theta=NULL)
{
	res = sapply(1:ncol(W), function(k) {
		if(!is.null(theta)) theta = lapply(theta, function(par) par[k])
			# get component-specific parameters to pass to estimation fn.
		thetahat.pvii(x, W[,k], aslist=FALSE, iter.max=0, theta=theta)
	})
	
	mean  = res["mean", ]
	scale = res["scale",]
	shape = res["shape",]
	
	if(aslist)
		namedList(mean, scale, shape)
	else
		namedVector(mean, scale, shape)
}
.thetahat_mvpvii <- function(X, W, aslist=TRUE, theta=NULL)
{
	kvec = 1:ncol(W)
	
	res = lapply(kvec, function(k) {
		if(!is.null(theta)) theta = lapply(theta, function(par)
				if(is.list(par)) par[[k]]
				else if(is.matrix(par)) par[k,]
				else par[k])
					# get component-specific parameters to pass on
		thetahat.mvpvii(X, W[,k], aslist=TRUE, iter.max=0, theta=theta)
	})
		
	mean  = t(sapply(kvec, function(k) res[[k]]$mean))
	scale = lapply(kvec, function(k) res[[k]]$scale)
	shape = sapply(kvec, function(k) res[[k]]$shape)
		# see .thetahat_mvnorm() above to understand the formation of these data structures
	
	if(aslist)
		namedList(mean, scale, shape)
	else
		namedVector(mean, scale, shape)
}
.thetahat_weisd <- function(x, W, aslist=TRUE, theta=NULL)
{
	res = sapply(1:ncol(W), function(k)
		thetahat.weisd(x, W[,k], aslist=FALSE, shape.guess=theta$shape[k]))
	
	sigma = res["shape",]
	delta = res["decay",]
		# the kth element of sigma is the shape for the kth component, and similarly for delta (decay)
	
	if(aslist)
		list(shape=sigma, decay=delta)
	else
		c(shape=sigma, decay=delta)
}
.thetahat_mvweisd <- function(X, W, aslist=TRUE, theta=NULL)
{
	kvec = 1:ncol(W)
	
	res = sapply(kvec,
	             function(k)
	             	thetahat.mvweisd(X, W[,k], aslist=TRUE, 
	             		shape.guess=theta$shape[k,]),
		         simplify=FALSE)
	
	sigma = t(sapply(kvec, function(k) res[[k]]$shape))
		# the kth row of sigma is the vector of the column shapes of X for the kth component
	delta = t(sapply(kvec, function(k) res[[k]]$decay))
		# the kth row of delta is the vector of the column decays of X for the kth component
	rho = lapply(kvec, function(k) res[[k]]$corr)
	
	list(shape=sigma, decay=delta, corr=rho)
}

# Build "Q" (transition) matrices under under local independence assumption (.qgen) or under linear dependence assumption (.qgenDep)
.qgen <- function(A, B)
{
	J = crossprod(A, B) # N f(A,B)
	
	J / rowSums(J) # f(B|A)
}
.qgenDep <- function(A, B)
{
	res = sapply(1:ncol(B), function(kb) coef(nnls(A, B[,kb]))) # N f(B|A)
	
	res / rowSums(res) # f(B|A)
}

# Build "V" (joint weight) matrices under local independence assumption (.vgen) or under linear dependence assumption (.vgenDep)
.vgen <- function(A, B) lapply(1:ncol(A), function(ka) A[,ka] * B)
.vgenDep <- function(A, B)
{
	R = .qgenDep(B, A) # f(A|B)
	tB = t(B)
	
	lapply(1:ncol(A), function(ka) t(tB * R[,ka]))
}

# Convert "V" (joint weight) matrix into MAP estimate.  TO DO:  this is horribly slow.
.v2map <- function(V)
{
	N = nrow(V[[1]])
	Ka = length(V)
	
	map.raw = lapply(1:N, function(n) { # MAP(V) indexed by n instead of ka
		res = t(sapply(1:Ka, function(ka) V[[ka]][n,]))
		1 * (res == max(res))
	})
	
	lapply(1:Ka, function(ka) t(sapply(1:N, function(N) map.raw[[N]][ka,])))
}