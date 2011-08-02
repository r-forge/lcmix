\name{logLik}

\alias{logLik}
\alias{logLik.mixmod}

\title{
Log-Likelihood for Mixture Models
}

\description{
S3 method to extract the log-likelihood for an object of class \code{mixmod} or \code{mdmixmod}.
}

\usage{
\method{logLik}{mixmod}(mod)
\method{logLik}{mdmixmod}(mod)
}

\arguments{
\item{mod}{an object of class \code{mixmod} or \code{mdmixmod}.}
}

\value{
Numeric.  Note that this is the log-likelihood for the parameters of the model given the observed data, \emph{not} the expected complete data log-likelihood as returned  by \code{\link{qval}}.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}}, \code{\link{reporting_functions}}, \code{\link{qfun}}.
}

\examples{\dontrun{		
data(exampleData1)
mod <- mixmod(exampleData1$X, 2)
logLik(mod) # [1] -1771.990
# logLik.mixmod(mod) is identical
}}

\keyword{models}
\keyword{cluster}
\keyword{classif}