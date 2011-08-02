\name{rank.sort.matrix}

\alias{rank.sort.matrix}

\title{
Order Matrix Rows
}

\description{
Puts a matrix in \emph{approximate} stochastic order, by row, using a function of the ranks in each row.
}

\usage{
rank.sort.matrix(X, decreasing=FALSE, fun=prod)
}

\arguments{
\item{X}{a numeric matrix.}
\item{decreasing}{logical; if \code{TRUE}, rows are given in descending order.}
\item{fun}{the function to apply to the column ranks of the matrix, either \code{prod} or \code{sum}.  Other functions may work, but produce bizarre results.}
}

\details{
Internally, \code{rank.sort.matrix} first ranks the values within the columns of \code{X}, then applies \code{fun} to the rows of the resulting rank matrix to determine the approximate stochastic order of the rows.
}

\value{
A matrix.
}

\author{
Daniel Dvorkin
}

\note{
The core algorithm in \code{rank.sort.matrix} is used to order components in the initialization procedure for \code{\link{mixmod}} and \code{\link{mdmixmod}}.
}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}}; \code{order}, \code{sort} in package \code{base}.
}

\examples{
set.seed(123)

X = matrix(rnorm(30), ncol=3)

rank.sort.matrix(X, TRUE) # approximately descending rows
#              [,1]       [,2]       [,3]
#  [1,]  1.71506499  1.7869131 -1.6866933
#  [2,]  0.46091621  0.4978505  0.8377870
#  [3,] -0.44566197 -0.4727914  1.2538149
#  [4,]  1.55870831  0.4007715 -1.0260044
#  [5,] -0.56047565  1.2240818 -1.0678237
#  [6,] -0.23017749  0.3598138 -0.2179749
#  [7,]  0.12928774 -0.5558411 -0.6250393
#  [8,]  0.07050839  0.1106827 -0.7288912
#  [9,] -0.68685285  0.7013559 -1.1381369
# [10,] -1.26506123 -1.9666172  0.1533731
}

\keyword{multivariate}
