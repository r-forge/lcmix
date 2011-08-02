
\name{reporting functions}
\alias{reporting_functions}
\alias{reporting functions}
\alias{bic}
\alias{bic.mixmod}
\alias{qval}
\alias{qval.mixmod}

\title{
Mixture Model Statistics
}

\description{
See Details for a description of the individual functions.
}

\usage{
bic(mod)
qval(mod)
}

\arguments{
\item{mod}{an object of class \code{\link{mixmod}} or \code{\link{mdmixmod}}.}
}

\details{
The function \code{bic} calculates the Bayes Information Criterion (Bayesian Information Criterion, BIC) measure for mixture models of class \code{mixmod} and \code{mdmixmod}.  BIC is defined here as \eqn{2 \mathcal{L} - |\Theta| \log N}{2*L - T*log(N)}, where \eqn{\mathcal{L}}{L} is the log-likelihood, \eqn{|\Theta|}{T} is the size of the parameter space, and \eqn{N}{N} is the size of the data.  Thus the model which provides the highest BIC should be selected.

The function \code{qval} provides the value of the Q-function, that is, the expected complete data (observed data and hidden components) log-likelihood for the parameters of the model.
}

\references{
McLachlan, G.J. and Thriyambakam, K.  (2008)  \emph{The EM Algorithm and Extensions}, John Wiley & Sons.
    
Schwarz, G. (1978)  Estimating the dimension of a model, \emph{The Annals of Statistics} \bold{6(2)}, 461--464.
}

\value{
Numeric.
}

\author{
Daniel Dvorkin
}

\note{
Some authors define BIC as the negative of the quantity given in Details. 
}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}}, \code{\link{qfun}}, \code{\link{logLik.mixmod}}.
}

\examples{\dontrun{
data(exampleData1)
mod <- mixmod(exampleData1$X, 2)
bic(mod) # -3578.518
qval(mod) # -1828.572
}}

\keyword{models}
\keyword{cluster}
\keyword{classif}