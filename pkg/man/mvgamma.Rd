\name{mvgamma}

\alias{mvgamma}
\alias{dmvgamma}
\alias{rmvgamma}

\title{
The Multivariate Gamma Distribution
}

\description{
Density and random generation functions for the multivariate gamma distribution constructed using a normal (Gaussian) copula.
}

\usage{
dmvgamma(X, shape, rate, corr=diag(ncol(X)), log=FALSE)
rmvgamma(N, shape=1, rate=1, corr=diag(length(shape)))
}

\arguments{
\item{X}{a numeric matrix of which each row represents an observation.}
\item{shape}{a vector of shape parameters for the marginal distributions of the columns of \code{X}.  If \code{length(shape) < ncol(X)}, the elements of \code{shape} will be recycled.  If \code{length(shape) > ncol(X)}, the \code{shape} vector will be truncated and a warning given.}
\item{rate}{a vector of rate parameters for the marginal distributions of the columns of \code{X}.  If \code{length(rate) < ncol(X)}, the elements of \code{rate} will be recycled.  If \code{length(rate) > ncol(X)}, the \code{rate} vector will be truncated and a warning given.}
\item{corr}{the correlation matrix.  See Details.}
\item{log}{logical; if \code{TRUE}, density is given as the log-density.}
\item{N}{number of vectors to simulate.}
}

\details{
The construction of multivariate distributions from univariate marginal distributions using normal copulas is discussed in Song (2000).
Briefly, given univariate marginal densities and the corresponding distribution functions (here, the gamma distribution), the standard normal quantiles of the values of the distribution
functions follow a multivariate standard normal distribution, that is, a multivariate normal distribution with marginal means of 0 and marginal variances of 1.
Thus the covariance matrix is referred to as the correlation matrix in this context.
} % line breaks inside paragraphs prevent Rd2txt from truncating lines

\value{
For \code{dmvgamma}, a vector of densities.  For \code{rmvgamma}, a vector with \code{N} rows and \code{ncol(corr)} columns representing a sample from the multivariate gamma distribution with the specified parameters.
}

\references{
Song, P.  (2000)  Multivariate dispersion models generated from Gaussian copula.  \emph{Scandinavian Journal of Statistics} \bold{27}, 305--320.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{mvnorm}}, \code{\link{mvweisd}} for related distributions; \code{\link{thetahat}} for parameter estimation.
}

\examples{
set.seed(123)

s <- 1:2
r <- 2:3
rho <- matrix(c(1, 0.5, 0.5, 1), ncol=2)

X <- rmvgamma(5, s, r, rho)
print(X)
#           [,1]      [,2]
# [1,] 0.1600585 0.4178748
# [2,] 1.3762076 0.7846068
# [3,] 0.6280634 1.6192241
# [4,] 0.3958020 0.2131819
# [5,] 0.1229824 0.3399217

dmvgamma(X, s, r, rho)
# [1] 1.92014830 0.10438760 0.08317813 0.79129064 2.32215333
}

\keyword{distribution}
\keyword{multivariate}
