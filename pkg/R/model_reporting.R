### model_reporting.R:  model reporting functions for the lcmix package.  A number of functions here are set as methods using the function setMethodS3(), which requires the package "R.methodsS3".  The package "ROCR" is also recommended, and is required for all functions under the header "ROC FUNCTIONS".

### LIKELIHOOD, ETC. FUNCTIONS
# TO DO:  there are a bunch of BIC, ABIC, deviance, log-likelihood, mutual information, etc. type model reporting functions sitting around in flydev_all_chrom/{include.R, mi.R, lcmix.OLD.R}.  What we need to do is decide on a final set of performance measures and include those here, without all the cruft.  2011-07-09 DD.

# BIC (Bayes information criterion) for an object of class "mixmod" or "mdmixmod".
setMethodS3("bic", "mixmod",
	function(mod) mod$stats[["bic"]],
	appendVarArgs=FALSE)

# Log-likelihood for an object of class "mixmod" or "mdmixmod".
setMethodS3("logLik", "mixmod",
	function(mod) mod$stats[["llik"]],
	appendVarArgs=FALSE)

# Q-function and Q-value for an object of class "mixmod" or "mdmixmod".
setMethodS3("qfun", "mixmod", function(mod)
{
	switch(class(mod)[1],
		   mixmod = .qfun_mixmod(mod$weights, mod$params, mod$pdfs),
		   mdmixmod = .qfun_mdmixmod(mod$weights, mod$params, mod$pdfs))
}, appendVarArgs=FALSE)
setMethodS3("qval", "mixmod",
	function(mod) sum(qfun(mod)),
	appendVarArgs=FALSE)

### PRINTING FUNCTIONS

# Print a representation of an object of class "mdmixmod".
setMethodS3("print", "mdmixmod", function(mod, ...)
{
	headstrbase = "%s %s mixture model %s\n"
	family = collapse.vec(sapply(mod$family, .pretty_family))
	distn = collapse.vec(mod$distn, quote=TRUE)
	headstr = upper.first(sprintf(headstrbase, mod$topology, family, distn))
	cat(headstr)
	
	parstrbase = "Data '%s' of size %i-by-%s fitted to %i %s components\n"
	D = collapse.vec(mod$D, space=FALSE)
	K = collapse.vec(mod$K, space=FALSE)
	parstr = sprintf(parstrbase, mod$dname, mod$N, D, mod$K0, K)
	cat(parstr)
	
	cat("Model statistics:\n")
	print(mod$stats, ...)
}, appendVarArgs=FALSE)

# Print a representation of an object of class "mixmod"
setMethodS3("print", "mixmod", function(mod, ...)
{
	family = .pretty_family(mod$family)
	
	headstrbase = "%s mixture model ('%s')\n"
	headstr = upper.first(sprintf(headstrbase, family, mod$distn))
	cat(headstr)
	
	parstrbase = "Data '%s' of %s fitted to %i components\n"
	moddim = ifelse(mod$D > 1,
		sprintf("size %i-by-%i", mod$N, mod$D),
		sprintf("length %i", mod$N))
	parstr = sprintf(parstrbase, mod$dname, moddim, mod$K)
	cat(parstr)
	
	cat("Model statistics:\n")
	print(mod$stats, ...)
}, appendVarArgs=FALSE)

### PLOTTING FUNCTIONS

# Plot the convergence of the log-likelihood for an object of class "mixmod" or "mdmixmod".  If show.qval is TRUE, also show Q-value.
convergence.plot <- function(mod, show.qval=FALSE, main="")
{
	if(show.qval) {
		ylim = range(mod$iteration.stats$llik, mod$iteration.stats$qval)
		ylab = "log-likelihood, Q-value"
	} else {
		ylim = range(mod$iteration.stats$llik)
		ylab = "log-likelihood"
	}
	
	plot(mod$iteration.stats$llik, type="l", xlab="iteration", ylab=ylab,
		ylim=ylim, main=main)
	
	if(show.qval)
		lines(mod$iteration.stats$qval, col="red")
}

### SUMMARY FUNCTIONS
# TO DO: see summary.mixmod() and print.summary.mixmod() in lcmix.OLD.R and decide if there's anything there worth keeping.  If you do decide to put them in, make sure they're registered as proper summary() and print() methods.  2011-07-09 DD.

### ROC FUNCTIONS

# Given one or more predictions (as a list of vectors, all of the same length N) and a set of labels (as a vector of length N, or a matrix or data frame with N rows and length(predictions) columns, or a list of vectors of the same dimension as predictions), fit ROC or quasi-ROC curves and calculate performance measures.  Return an object of class "multiROC", a list with the following elements:
# $pred:  a list of objects of class "prediction"
# $perf:  a list of objects of class "performance"
# $meas:  a list of vectors with named elements "auc" and "dca" (ROC performance measures)
# $quasi:  Boolean, TRUE if this is for quasi-ROC, FALSE otherwise
multiROC <- function(predictions, labels, quasi=FALSE)
{
	require("ROCR")
	
	if(is.matrix(predictions)) predictions = as.data.frame(predictions)
	
	if(is.null(names(predictions)))
		names(predictions) = sprintf("P%i", 1:length(predictions))
	
	if(is.matrix(labels)) labels = as.data.frame(labels)
	
	if(is.list(labels)) {
		names(labels) = names(predictions)
		pred = lapply(names(predictions), function(pname)
				prediction(predictions[[pname]], labels[[pname]]))
		names(pred) = names(predictions)
	} else {
		pred = lapply(predictions, function(p) prediction(p, labels))
	}
	
	measure = "tpr"
	x.measure = ifelse(quasi, "rpp", "fpr")
	perf = lapply(pred, performance, measure=measure, x.measure=x.measure)
	
	meas = lapply(perf,
		function(p)
		{
			x = p@x.values[[1]]
			y = p@y.values[[1]]
			auc = sum(y * diff(c(x, 1)))
			dca = sqrt(min(x^2 + (y-1)^2))
			c(auc=auc, dca=dca)
		})
	
	retn = list(pred=pred, perf=perf, meas=meas, quasi=quasi)
	class(retn) = "multiROC"
	return(retn)
}

# Plot an object of class "multiROC".
setMethodS3("plot", "multiROC", function(mroc,
	do.legend=TRUE, do.auc=TRUE, do.dca=TRUE,
	xlab=ifelse(mroc$quasi, "all positive rate", "false positive rate"),
	ylab="true positive rate", main="", colorplus=0, colors=colorplus+(1:length(mroc$perf)), ...)
{
	require("ROCR")
	
	pnames = names(mroc$perf)
	nnames = length(pnames)
	
	plot(mroc$perf[[pnames[1]]], xlab=xlab, ylab=ylab, main=main, col=colors[1], 
		...)
	
	if(length(pnames) > 1)
		for(i in 2:length(pnames))
			plot(mroc$perf[[pnames[i]]], col=colors[i], add=TRUE, ...)
	
	if(do.legend) {
		ltext = pnames
		if(do.auc) {
			auc = sapply(mroc$meas, function(m) m[["auc"]])
			ltext = paste(ltext, sprintf(", AUC=%.3f", auc), sep="")
		}
		if(do.dca) {
			dca = sapply(mroc$meas, function(m) m[["dca"]])
			ltext = paste(ltext, sprintf(", DCA=%.3f", dca), sep="")
		}
		legend("bottomright", lty=1, col=colors, legend=ltext)
	}
}, appendVarArgs=FALSE)

# Plot an object of class "multiROC" in black and white.  TO DO:  it would be good to eventually just make BW plotting an option in plot.multiROC.  2011-07-09 DD.
plot.multiROC.bw <- function(mroc, do.legend=TRUE, do.auc=TRUE, do.dca=TRUE,
	xlab=ifelse(mroc$quasi, "all positive rate", "false positive rate"),
	ylab="true positive rate", main="", ltys=1:length(mroc$perf), ...)
{
	require("ROCR")
	
	pnames = names(mroc$perf)
	nnames = length(pnames)
	
	plot(mroc$perf[[pnames[1]]], xlab=xlab, ylab=ylab, main=main, lty=ltys[1], 
		...)
	
	if(length(pnames) > 1)
		for(i in 2:length(pnames))
			plot(mroc$perf[[pnames[i]]], lty=ltys[i],, lwd=1.5, add=TRUE, ...)
	
	if(do.legend) {
		ltext = pnames
		if(do.auc) {
			auc = sapply(mroc$meas, function(m) m[["auc"]])
			ltext = paste(ltext, sprintf(", AUC=%.3f", auc), sep="")
		}
		if(do.dca) {
			dca = sapply(mroc$meas, function(m) m[["dca"]])
			ltext = paste(ltext, sprintf(", DCA=%.3f", dca), sep="")
		}
		legend("bottomright", lty=ltys, legend=ltext)
	}
}
