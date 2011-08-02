\name{convergence.plot}

\alias{convergence.plot}

\concept{convergence}

\title{
Show EM Algorithm Convergence
}

\description{
Plot convergence of the log-likelihood to a maximum during execution of the EM algorithm for fitting mixture models.  Optionally also show the Q-value, that is, the expected complete data log-likelihood.
}
\usage{
convergence.plot(mod, show.qval=FALSE, main="")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{mod}{an object of class \code{mixmod} or \code{mdmixmod}.}
\item{show.qval}{logical; if \code{TRUE}, the Q-value will be plotted as well as the log-likelihood.}
\item{main}{plot title.}
}

\details{
Log-likelihood is plotted vs. iteration in black.  If \code{show.qval} is \code{TRUE}, Q-value is plotted in red.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}}, \code{\link{exampleData}}, \code{\link{multiROC}}; \code{plot} in package \code{graphics}.
}

\examples{\dontrun{
data(exampleData1)
mod1 = mixmod(exampleData1$X, 2)
convergence.plot(mod1, TRUE, main="log-likelihood and Q-value")
}}

\keyword{hplot}
