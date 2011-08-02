\name{thetahat}

\alias{thetahat}
\alias{thetahat.norm}
\alias{thetahat.mvnorm}
\alias{thetahat.weisd}
\alias{thetahat.mvweisd}
\alias{thetahat.gamma}
\alias{thetahat.mvgamma}

\concept{parameter estimation}
\concept{maximum likelihood estimation}
\concept{maximum likelihood estimator}

\title{
Parameter Estimation for Common Distributions
}

\description{
Functions to estimate the parameters of distributions used in the \code{lcmix} package by maximum likelihood estimation.
}

\usage{
thetahat.norm(x, w=1, aslist=TRUE)
thetahat.mvnorm(X, w=1, aslist=TRUE)
thetahat.weisd(x, w=1, aslist=TRUE)
thetahat.mvweisd(X, w=1, aslist=TRUE)
thetahat.gamma(x, w=1, aslist=TRUE)
thetahat.mvgamma(X, w=1, aslist=TRUE)
}

\arguments{
\item{x}{a numeric vector, with each element representing an observation.}
\item{X}{a numeric matrix or data frame, with each row representing an observation.}
\item{w}{a numeric vector of non-negative values; weights to be applied to the elements of \code{x} or rows of \code{X} in parameter estimation.}
\item{aslist}{logical; if \code{TRUE}, estimated parameters will be returned as a named list.}
}

\details{
By default, estimation is unweighted.  If weights are given, the length of the weights should generally be the same as the length of the data, that is, \code{length(w) = length(x)} or \code{length(w) = nrow(X)}.  If weights are shorter than the data, they will be recycled.  If longer than the data, weights will be truncated and a warning will be given.

If \code{aslist} is \code{FALSE}, parameters will be returned as a vector with named elements if possible.  However, multivariate distributions (\code{mv*}) often have matrix parameters and using \code{aslist=FALSE} for these distributions may produce results which are difficult to interpret.
}

\value{
A named list, a vector, or a matrix, depending on the distribution and options used, of estimated parameters.
}

\references{
Lawless, J.F.  (1982)  \emph{Statistical Models and Methods for Lifetime Data}, John Wiley & Sons.

Joe, H. and Xu, J.  (1996)  The estimation method of inference functions for margins for multivariate models.  Technical Report 166, University of British Columbia, Department of Statistics.
}

\author{
Daniel Dvorkin
}

\note{
For distributions without closed-form MLEs, a method-of-moments approach is used to obtain initial parameter estimates and then numerical optimization of the likelihood function is performed using \code{nlminb} in package \code{stats}.  See Lawless (1982) for issues in parameter estimation for the Weibull and gamma distributions.

For multivariate distributions other than the normal, the correlation matrix for the copula is estimated using the inference function for margins procedure of Joe and Xu (1996).
}

\seealso{
\code{\link{mvnorm}}, \code{\link{weisd}}, \code{\link{mvweisd}}, \code{\link{mvgamma}}; \code{nlminb}, \code{Normal}, \code{Weibull}, \code{Gamma} in package \code{stats}; \code{fitdistr} in package \code{MASS} which performs parameter estimation using different numerical optimization functions.
}

\examples{
set.seed(123)
x <- rweisd(1000, 2, 3)
thetahat.weisd(x, aslist=FALSE)
#    shape    decay
# 2.015511 2.996586
}

\keyword{distribution}
\keyword{univar}
\keyword{multivariate}