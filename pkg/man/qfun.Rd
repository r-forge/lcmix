\name{qfun}

\alias{qfun}
\alias{qfun.mixmod}

\title{
Q-Function Values for Mixture Models
}

\description{
Calculate hidden and observed data components of the Q-function (expected complete data log-likelihood).
}

\usage{
qfun(mod)
}

\arguments{
\item{mod}{an object of class \code{\link{mixmod}} or \code{\link{mdmixmod}}.}
}

\details{
The Q-function is the expected value of the complete data (observed data and hidden components) log-likelihood for the parameters of the model.
The EM algorithm used in model fitting attempts to maximize this function.
} % line breaks inside paragraphs prevent Rd2txt from truncating lines

\value{
A vector with the following named components:
\item{hidden}{the hidden component log-likelihood.}
\item{observed}{the observed data log-likelihood.}
}

\references{
McLachlan, G.J. and Thriyambakam, K.  (2008)  \emph{The EM Algorithm and Extensions}, John Wiley & Sons.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}}, \code{\link{reporting_functions}}, \code{\link{logLik.mixmod}}.
}

\examples{\dontrun{		
data(exampleData1)
mod <- mixmod(exampleData1$X, 2)
qfun(mod)
#    hidden  observed 
# -308.2579 -1520.3143
}}

\keyword{models}
\keyword{cluster}
\keyword{classif}