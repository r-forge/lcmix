### simulation.R:  simulation functions for the lcmix package.

# Simulate multiple mixture data of the given size, distributions, parameters, and topology; params should be as the $params element of the return value of mdmixmod().  Return a list with elements $X (a list of observed data), $Y (a list of hidden components corresponding to the observed data), and $Y0 (top-level hidden components), with attribute topology (a string giving the topology, i.e. "layered" or "chained".)
simulate.mdmixdata <- function(N, distn, params, topology=LC_TOPOLOGY)
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
			mvgamma = ncol(params$observed[[z]]$shape)))
	
	K0 = length(params$hidden$prob0)
	K = sapply(params$hidden$cprob, ncol)
	
	## simulate top-level hidden component
	
	Y0 = rcategorical(N, params$hidden$prob0)
	
	## simulate lower-level hidden components
	
	if(topology == "layered") {
	
		Y = lapply(zvec, function(z)
		{
			res = vector(mode="numeric", length=N)
			for(k0 in 1:K0)
			{
				isk0 = (Y0 == k0)
				nisk0 = sum(isk0)
				res[isk0] = rcategorical(nisk0, params$hidden$cprob[[z]][k0,])
			}
			return(res)
		})
	
	} else if(topology == "chained") {
		
		Y = vector(mode="list", length=Z)
		
		Y[[1]] = vector(mode="numeric", length=N)
		for(k0 in 1:K0)
		{
			isk0 = (Y0 == k0)
			nisk0 = sum(isk0)
			Y[[1]][isk0] = rcategorical(nisk0, params$hidden$cprob[[1]][k0,])
		}
		
		if(Z > 1) for(z in 2:Z)
		{
			Y[[z]] = vector(mode="numeric", length=N)
			for(kzm1 in 1:K[[z-1]])
			{
				iskzm1 = (Y[[z-1]] == kzm1)
				niskzm1 = sum(iskzm1)
				Y[[z]][iskzm1] = rcategorical(niskzm1, 
					params$hidden$cprob[[z]][kzm1,])
			}
		}
	
		names(Y) = zvec
	
	} # "else":  unsupported topology, to have been caught by "match.arg" above
	
	## simulate observed data
	
	X = lapply(zvec, function(z)
	{
		par = params$observed[[z]]
		
		if(D[[z]] == 1) {
			res = vector(mode="numeric", length=N)
		} else {
			res = matrix(0, nrow=N, ncol=D[[z]])
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
				gamma = rgamma(niskz, par$shape[kz], par$rate[kz]))
			else res[iskz,] = switch(distn[[z]],
				mvnorm = rmvnorm(niskz, par$mean[kz,], par$cov[[kz]]),
				mvweisd = rmvweisd(niskz, par$shape[kz,], par$decay[kz,], 
					par$corr[[kz]]),
				mvgamma = rmvgamma(niskz, par$shape[kz,], par$rate[kz,], 
					par$corr[[kz]]))
		}
		
		return(res)
	})
	
	## package and return
	
	retn = list(X=X, Y=Y, Y0=Y0)
	attr(retn, "topology") = topology
	return(retn)
}
			
# Simulate mixture data of the given data size, distribution, and parameters; params should be as the $params element of the return value of mixmod(). Return a list with elements $X (observed data) and $Y (hidden component).
simulate.mixdata <- function(N, distn, params)
{
	## simulate hidden component
	
	Y = rcategorical(N, params$hidden$prob)
	
	## simulate observed data
	
	# gather statistics and set up observed data
	
	D = switch(distn,
			norm    = 1,
			mvnorm  = ncol(params$observed$mean),
			weisd   = 1,
			mvweisd = ncol(params$observed$shape),
			gamma   = 1,
			mvgamma = ncol(params$observed$shape))
	
	if(D == 1) {
		X = vector(mode="numeric", length=N)
	} else {
		X = matrix(0, nrow=N, ncol=D)
		if(distn == "mvnorm")
			colnames(X) = colnames(params$observed$mean)
		else if(distn %in% c("mvweisd", "mvgamma"))
			colnames(X) = colnames(params$observed$shape)
	}
	
	## simulate observed data
	
	par = params$observed
	
	for(k in 1:length(params$hidden$prob))
	{
		isk = (Y == k)
		nisk = sum(isk)
		
		if(D == 1) X[isk] = { # TO DO:  replace these chains of "if-else-if" with "switch" statements a la simulate.mdmixdata().  2011-07-17 DD.
			if(distn == "norm")
				rnorm(nisk, par$mean[k], par$sd[k])
			else if(distn == "weisd")
				rweisd(nisk, par$shape[k], par$decay[k])
			else if(distn == "gamma")
				rgamma(nisk, par$shape[k], par$rate[k])
		} else X[isk,] = {
			if(distn == "mvnorm")
				rmvnorm(nisk, par$mean[k,], par$cov[[k]])
			else if(distn == "mvweisd")
				rmvweisd(nisk, par$shape[k,], par$decay[k,], par$corr[[k]])
			else if(distn == "mvgamma")
				rmvgamma(nisk, par$shape[k,], par$rate[k,], par$corr[[k]])
		}
	}
				
	## package and return
	
	list(X=X, Y=Y)
}
	
# simulate a model with the same size, distribution, and parameters as the fitted model "mod", an object of class "mdmixmod"
setMethodS3("simulate.from.fit", "mdmixmod", function(mod)
	simulate.mdmixdata(mod$N, mod$distn, mod$params, mod$topology),
	appendVarArgs=FALSE)

	
# simulate a model with the same size, distribution, and parameters as the fitted model "mod", an object of class "mixmod"
setMethodS3("simulate.from.fit", "mixmod", function(mod)
	simulate.mixdata(mod$N, mod$distn, mod$params),
	appendVarArgs=FALSE)
