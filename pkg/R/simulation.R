### simulation.R:  simulation functions for the lcmix package.

# resample from data to create a data set with the same size, distribution, and hidden parameters as the fitted model "x", an object of class "mixmod" or "mdmixmod"
setMethodS3("resampleFromFit", "mixmod",
function(x, n=x$N, replace=TRUE, hidpar=x$params$hidden, ...)
{
	# simulate hidden data
	Y = .simulate_hidden_mixdata(n, hidpar)
	
	# simulate observed data
	x$X = as.matrix(x$X)
	X = matrix(nrow=n, ncol=ncol(x$X))
	for(k in x$kvec)
	{
		isk = (Y == k)
		nisk = sum(isk)
		prob = x$weights$W[,k]
		X[isk,] = rowSample(x$X, nisk, replace, prob)
	}
	X = drop(X)
	
	# package it up, send it back
	namedList(X, Y)
})
setMethodS3("resampleFromFit", "mdmixmod",
function(x, n=x$N, replace=TRUE, hidpar=x$params$hidden,
	topology=x$topology, ...)
{
	# simulate hidden data
	hidden = .simulate_hidden_mdmixdata(n, hidpar, topology)
	Y = hidden$Y
	Y0 = hidden$Y0
	
	# simulate observed data
	X = lapply(x$zvec, function(z)
	{
		xXz = as.matrix(x$X[[z]])
		Xz = matrix(nrow=n, ncol=ncol(xXz))
		for(kz in x$kvec[[z]])
		{
			iskz = (Y[[z]] == kz)
			niskz = sum(iskz)
			prob = x$weights$W[[z]][,kz]
			Xz[iskz,] = rowSample(xXz, niskz, replace, prob)
		}
		drop(Xz)
	})
	
	# package it up, send it back
	namedList(X, Y, Y0)
})

# simulate data with similar performance to the Ci data
simulateCiTypeData <- function(n=10000, par=LC_SIMPAR, topology=LC_TOPOLOGY)
{
	topology = match.arg(topology)
	
	if(topology == "layered") {
		par$hidden$cprob = lapply(par$hidden$probz, function(p)
			.dependentDiscreteCpdf(par$hidden$prob0, p, dependence=0.60))
	} else { # topology == "chained"
		par$hidden$cprob = list(
			binding = .dependentDiscreteCpdf(par$hidden$prob0,
				par$hidden$probz$binding, dependence=0.85),
			expression = .dependentDiscreteCpdf(par$hidden$probz$binding,
				par$hidden$probz$expression, dependence=0.75),
			conservation = .dependentDiscreteCpdf(par$hidden$probz$expression,
				par$hidden$probz$conservation, dependence=0.75))
	}
	
	simulateMdmixdata(n, c("norm", "mvnorm", "norm"), par, topology)
}

# simulate a data set with the same size, distribution, and parameters as the fitted model "x", an object of class "mixmod" or "mdmixmod"
setMethodS3("simulateFromFit", "mixmod",
function(x, n=x$N, ...)
	simulateMixdata(n, x$distn, x$params))
setMethodS3("simulateFromFit", "mdmixmod",
function(x, n=x$N, ...)
	simulateMdmixdata(n, x$distn, x$params, x$topology))

	
# Simulate multiple mixture data of the given size, distributions, parameters, and topology; params should be as the $params element of the return value of mdmixmod().  Return a list with elements $X (a list of observed data), $Y (a list of hidden components corresponding to the observed data), and $Y0 (top-level hidden components), with attributes @topology (a string giving the topology, i.e. "layered" or "chained") and @params (the params used in the simulation.)
simulateMdmixdata <- function(n, distn, params, topology=LC_TOPOLOGY)
{
	## gather statistics, munge input PRN
		
	topology = match.arg(topology)
	
	Z = length(params$hidden$cprob)
	zvec = names(params$hidden$cprob)
	if(is.null(zvec))
		zvec = paste("X", 1:Z, sep="")
	names(params$hidden$cprob) = names(params$observed) = names(distn) = zvec
		# just to make sure everything is using the same names
	names(zvec) = zvec
		# so results of lapply-type functions have proper names
	
	D = sapply(zvec, function(z) switch(distn[[z]],
			norm    = 1,
			mvnorm  = ncol(params$observed[[z]]$mean),
			weisd   = 1,
			mvweisd = ncol(params$observed[[z]]$shape),
			gamma   = 1,
			mvgamma = ncol(params$observed[[z]]$shape),
			pvii    = 1,
			mvpvii = ncol(params$observed[[z]]$mean)))
			
	
	K0 = length(params$hidden$prob0)
	K = sapply(params$hidden$cprob, ncol)
	
	## simulate hidden data
	
	hidden = .simulate_hidden_mdmixdata(n, params$hidden, topology)
	Y = hidden$Y
	Y0 = hidden$Y0
	
	## simulate observed data
	
	X = lapply(zvec, function(z)
	{
		par = params$observed[[z]]
		
		if(D[[z]] == 1) {
			res = vector(mode="numeric", length=n)
		} else {
			res = matrix(0, nrow=n, ncol=D[[z]])
			if(distn[[z]] == "mvnorm")
				colnames(res) = colnames(params$observed$mean)
			else if(distn[[z]] %in% c("mvweisd", "mvgamma"))
				colnames(res) = colnames(params$observed$shape)
		}
		
		for(kz in 1:K[[z]])
		{
			iskz = (Y[[z]] == kz)
			niskz = sum(iskz)
			
			if(D[[z]] == 1) res[iskz] = switch(distn[[z]],
				norm = rnorm(niskz, par$mean[kz], par$sd[kz]),
				weisd = rweisd(niskz, par$shape[kz], par$decay[kz]),
				gamma = rgamma(niskz, par$shape[kz], par$rate[kz]),
				pvii = rpvii(niskz, par$mean[kz], par$scale[kz], 
					par$shape[kz]),
				stop(sprintf("unknown distribution '%s'", distn)))
			else res[iskz,] = switch(distn[[z]],
				mvnorm = rmvnorm(niskz, par$mean[kz,], par$cov[[kz]]),
				mvweisd = rmvweisd(niskz, par$shape[kz,], par$decay[kz,], 
					par$corr[[kz]]),
				mvgamma = rmvgamma(niskz, par$shape[kz,], par$rate[kz,], 
					par$corr[[kz]]), # MVPVII TO DO
				stop(sprintf("unknown distribution '%s'", distn)))
		}
		
		return(res)
	})
	
	## package and return
	
	retn = list(X=X, Y=Y, Y0=Y0)
	attr(retn, "n") = n
	attr(retn, "distn") = distn
	attr(retn, "params") = params
	attr(retn, "topology") = topology
	
	return(retn)
}
			
# Simulate mixture data of the given data size, distribution, and parameters; params should be as the $params element of the return value of mixmod(). Return a list with elements $X (observed data) and $Y (hidden component) and attribute @params (the parameters used in the simulation).
simulateMixdata <- function(n, distn, params)
{
	## simulate hidden data
	
	Y = .simulate_hidden_mixdata(n, params$hidden)
	
	## simulate observed data
	
	# gather statistics and set up observed data
	
	D = switch(distn,
			norm    = 1,
			mvnorm  = ncol(params$observed$mean),
			weisd   = 1,
			mvweisd = ncol(params$observed$shape),
			gamma   = 1,
			mvgamma = ncol(params$observed$shape),
			pvii    = 1,
			mvpvii  = ncol(params$observed$mean),
			stop(sprintf("distn '%s' not yet supported", distn)))
	
	if(D == 1) {
		X = vector(mode="numeric", length=n)
	} else {
		X = matrix(0, nrow=n, ncol=D)
		if(distn %in% c("mvnorm", "mvpvii"))
			colnames(X) = colnames(params$observed$mean)
		else if(distn %in% c("mvweisd", "mvgamma"))
			colnames(X) = colnames(params$observed$shape)
	}
	
	# simulate
	
	par = params$observed
	
	for(k in 1:length(params$hidden$prob))
	{
		isk = (Y == k)
		nisk = sum(isk)
		
		if(D == 1) X[isk] = { # TO DO:  replace these chains of "if-else-if" with "switch" statements a la simulateMdmixdata().
			if(distn == "norm")
				rnorm(nisk, par$mean[k], par$sd[k])
			else if(distn == "weisd")
				rweisd(nisk, par$shape[k], par$decay[k])
			else if(distn == "gamma")
				rgamma(nisk, par$shape[k], par$rate[k])
			else if(distn == "pvii")
				rpvii(nisk, par$mean[k], par$scale[k], par$shape[k])
		} else X[isk,] = {
			if(distn == "mvnorm")
				rmvnorm(nisk, par$mean[k,], par$cov[[k]])
			else if(distn == "mvweisd")
				rmvweisd(nisk, par$shape[k,], par$decay[k,], par$corr[[k]])
			else if(distn == "mvgamma")
				rmvgamma(nisk, par$shape[k,], par$rate[k,], par$corr[[k]])
			else if(distn == "mvpvii")
				rmvpvii(nisk, par$mean[k,], par$scale[[k]], par$shape[k])
		}
	}
				
	## package and return
	
	retn = list(X=X, Y=Y)
	attr(retn, "n") = n
	attr(retn, "distn") = distn
	attr(retn, "params") = params
	return(retn)
}