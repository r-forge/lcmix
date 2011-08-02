\name{exampleData}

\alias{exampleData}
\alias{exampleData1}
\alias{exampleData2}
\alias{exampleData3}
\alias{exampleData4}

\docType{data}

\title{
Example Data Sets
}

\description{
Examples of the types of data to which \code{\link{mixmod}} and \code{\link{mdmixmod}} may be used to fit mixture models, with complete (observed and hidden) data.
}

\format{
For \code{exampleData1} and \code{exampleData2}, a list with two elements, \code{X} and \code{Y}.  In \code{exampleData1}, \code{X} is a numeric vector, and \code{Y} is an integer vector of the same length as \code{X} containing values of either 1 or 2.  In \code{exampleData2}, \code{X} is a numeric matrix, and \code{Y} is an integer vector such that \code{length(Y) = nrow(X)}, containing values of 1, 2, or 3.

For \code{exampleData3} and \code{exampleData4}, a list with three elements, \code{X}, \code{Y}, and \code{Y0}, as described below:

\code{X} is a list with two elements, \code{ed.1} and \code{ed.2}, where \code{.} represents either \code{'3'} or \code{'4'}.  For example, the elements of \code{exampleData3$X} are \code{ed31} and \code{ed32}.  The first element, \code{ed.1}, is a numeric vector; the second, \code{ed.2}, is a numeric matrix such that \code{nrow(ed.2) = length(ed.1)}.

\code{Y} is a list with two elements, having the same names as the elements of \code{X}:  \code{ed.1} and \code{ed.2} are both integer vectors such that \code{length(Y$ed.1) = length(Y$ed.2) = length(X$ed.1)}.  The values in \code{ed.1} are 1 or 2, and the values in \code{ed.2} are 1, 2, or 3.

code{Y0} is an integer vector such that \code{length(Y0) = length(X$ed.1)}, containing values of either 1 or 2.

Both \code{exampleData3} and \code{exampleData4} also have the attribute \code{topology}, which is the string \code{"layered"} for \code{exampleData3}, and \code{"chained"} for \code{exampleData4}.
}

\details{
In all of the \code{exampleData.} sets (where \code{.} represents \code{'1'}, \code{'2'}, \code{'3'}, or \code{'4'}) the \code{X} element is the observed data, subject to model fitting with \code{\link{mixmod}} or \code{\link{mdmixmod}}, while the \code{Y} and \code{Y0} elements are the hidden components which actually generate the observed data, that is, specify the distribution from which the observed data are drawn.

\code{exampleData1} contains data generated from a 2-component univariate normal mixture distribution, while \code{exampleData2} contains data generated from a 3-component multivariate Weibull (shape-decay) mixture distribution.

\code{exampleData3} contains multiple data generated from a layered model in which the observed data come from a 2-component univariate normal mixture and a 3-component multivariate normal mixture, with the component distributions themselves being from a 2-component top-level mixture.  In other words, this data set represents a model in which the categorical random variable \eqn{\mathcal Y_0}{Y_0} generates the categorical random variables \eqn{\mathcal Y_1}{Y_1} and \eqn{\mathcal Y_2}{Y_2}, which in turn generate the univariate normal random variable \eqn{\mathcal X_1}{X_1} and the multivariate normal random variable \eqn{\mathcal X_2}{X_2}, respectively.

\code{exampleData4} contains multiple data generated from a chained model in which the observed data come from a 2-component univariate Weibull (shape-decay) mixture and a 3-component multivariate gamma mixture distribution, with the component distributions themselves being part of a chain of mixtures initiated by a top-level component.  In other words, this data set represents a model in which the categorical random variable \eqn{\mathcal Y_0}{Y_0} generates the categorical random variable \eqn{\mathcal Y_1}{Y_1}, which in turn generates the univariate Weibull random variable \eqn{\mathcal X_1}{X_1} and the categorical random variable \eqn{\mathcal Y_2}{Y_2}, and \eqn{\mathcal Y_2}{Y_2} generates the multivariate gamma random variable \eqn{\mathcal X_2}{X_2}.

For both \code{exampleData3} and \code{exampleData4}, the elements of \code{Y} give the components which specify the distributions from which the elements of \code{X} are drawn.  Also, in \code{exampleData3}, \code{Y0} gives the components which specify the distributions from which the elements of \code{Y} are drawn.  And in \code{exampleData4}, \code{Y0} gives the components which specify the distribution from which \code{Y$ed41} is drawn; in turn, \code{Y$ed41} gives the components which specify the distribution from which \code{Y$ed42} is drawn.

The \code{topology} attribute of \code{exampleData3} and \code{exampleData4} gives the model topology under which the data were generated. 
}

\examples{\dontrun{
data(exampleData1)
mod1 = mixmod(exampleData1$X, 2)

data(exampleData2)
mod2 = mixmod(exampleData2$X, 3, family="weibull")

data(exampleData3)
mod3 = mdmixmod(exampleData3$X, K=c(2,3), K0=2)

data(exampleData4)
mod4 = mdmixmod(exampleData4$X, K=c(2,3), K0=2,
                family=c("wei", "gam"), topology="chained")
    # note use of partial family names is acceptable
}}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}} for the functions used in Examples; \code{\link{simulation}} for the functions used to generate the example data sets; \code{\link{mvnorm}}, \code{\link{mvweisd}}, \code{\link{weisd}}, \code{\link{mvgamma}} for the relevant distributions; \code{\link{multiROC}} for performance evaluations making use of complete data.
}

\keyword{datasets}