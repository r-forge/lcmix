\name{euclidean}

\alias{euclidean}

\title{
Euclidean Distance
}

\description{
Returns the squared Euclidean distance between all rows in \code{X} and the vector \code{v}.
}

\usage{
euclidean(X, v, ...)
}

\arguments{
\item{X}{a matrix or data frame, or a vector to be interpreted as a one-row matrix.}
\item{v}{a vector such that \code{length(v) = ncol(X)}.}
\item{\dots}{no effect; included for signature compatibility with \code{mahalanobis}.}
}

\value{
A vector of which the \code{i}th element is the \emph{squared} Euclidean distance between \code{X[i,]} and \code{v}.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{mahalanobis} in package \code{stats}.
}

\examples{
set.seed(123)
Z <- matrix(rnorm(12), ncol=3)
euclidean(Z, rep(0,3))
# [1] 0.8026151 3.1930442 4.1403916 1.7348174
}

\keyword{array}
\keyword{arith}
