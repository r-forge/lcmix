### visible_utilities.R:  visible utility functions for the lcmix package

# Return a printable "collapsed" version of a vector.  If "space" is TRUE, then the elements will be set off by spaces.  If quote is TRUE, then single quotes will be added around the elements of the vector
collapse.vec <- function(x, space=TRUE, quote=FALSE)
{	
	collapse = ifelse(space, ", ", ",")
	if(quote) x = sprintf("'%s'", x)
	
	sprintf("(%s)", paste(x, collapse=collapse))
}

# Return an N-length vector containing the squared Euclidean distances between all rows of the N-by-D matrix X and the D-length vector v.  This is by analogy to the built-in function mahalanobis() which returns the squared Mahalanobis distance.  The "..." argument exists solely to allow euclidean() to be called with additional named arguments, such as the "cov" and "inverted" arguments to mahalanobis(); these arguments will have no effect on the value returned by euclidean().
euclidean <- function(X, v, ...)
{
	if(is.vector(X)) X = matrix(X, nrow=1)
	
	colSums((t(as.matrix(X)) - v)^2)
}

# replace all values in X less than min with min, and all values greater than max with max
force.range <- function(X, min=-Inf, max=Inf)
{
	X[X < min] = min
	X[X > max] = max
	
	return(X)
}

# Return the MLE of the covariance of X assuming that mu is the mean.  By default, if X is univariate, a scalar will be returned, but if simplify is FALSE, then a matrix will be returned regardless.
mlecov <- function(X, mean=colMeans(as.matrix(X)), simplify=TRUE)
{
	X = as.matrix(X)
	retn = .mleTcov(t(X), mean, nrow(X))
	if(simplify) drop(retn) else retn
}

# Return the MLE of the standard deviation of X, or if X is a matrix or data frame, the vector of column standard deviations.
mlesd <- function(X) drop(sqrt(diag(mlecov(X, simplify=F))))

# put a matrix in rough stochastic order
rank.sort.matrix <- function(X, decreasing=FALSE, fun=prod)
{
	if(decreasing)
		ranks = apply(-X, 2, rank)
	else
		ranks = apply(X, 2, rank)
	
	rowranks = apply(ranks, 1, fun)
	
	X[order(rowranks),]
}

# Standardize data; if X is a vector, return the standardized vector; if X is a matrix or data frame, return with standardized columns  In the case of method="sd", this means dividing data by its standard deviation; in the case of method="mean", this means dividing data by its mean.  If mle is TRUE, then the maximum likelihood estimator of standard deviation rather than the unbiased estimator will be used.  (The MLE argument has no effect if method="mean".)  if center is TRUE, the data will also be centered to have mean (or column means) 0.  Thus "standardize(X, 'sd', FALSE, TRUE)" is equivalent to "scale(X)".  The return value will have attribute standardized, which is a list with elements $scale and $center, which have the same meanings as scaled:scale and scaled:center in the return value attributes of scale().
standardize <- function(X, method=c("sd", "mean"), mle=TRUE, center=FALSE)
{
	method = match.arg(method)
	
	asmat = is.matrix(X) || is.data.frame(X)
	if(asmat) tX = t(X) # we will use this repeatedly below
	
	mu = if(asmat) colMeans(X) else mean(X)
	if(method == "sd") sigma = if(mle) mlesd(X) else sd(X)
	
	if(center) {
		center = mu
		if(asmat) {
			tX = tX-center
			X = t(tX)
		} else X = X-center
	} else center = NULL
	
	scale = switch(method, sd=sigma, mean=mu)
	X = if(asmat) t(tX/scale) else X/scale
	
	retn = X
	attr(retn, "standardized") = list(scale=scale, center=center)
	return(retn)
}
standardise <- standardize

# Return k evenly spaced quantiles of X (or quantiles of the columns of X, if X is a matrix or data frame) which are by default ordered from greatest to least.
uniform.quantiles <- function(X, K, decreasing=TRUE)
{
	if(decreasing)
		ssq = K:1
	else
		ssq = 1:K
	
	probs = (ssq - 0.5) / K
	
	if(is.null(dim(X))) # vector
		quantile(X, probs, names=FALSE)
	else # matrix or data frame
		apply(X, 2, quantile, probs=probs, names=FALSE)
}

# Given a string, change the first letter in the string (if any) to uppercase
upper.first <- function(s)
{
	pos = regexpr("[A-za-z]", s)
		# check for upper-case already there so as not to capitalize the next letter (if any) after the first upper-case one
	if(pos > 0)
		substr(s, pos, pos) = toupper(substr(s, pos, pos))
	
	return(s)
}