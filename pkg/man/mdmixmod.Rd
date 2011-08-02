\name{mdmixmod}

\alias{mdmixmod}
\alias{print.mdmixmod}

\concept{parameter estimation}
\concept{maximum likelihood estimation}
\concept{maximum likelihood estimator}
\concept{mixture model estimation}
\concept{EM algorithm}
\concept{expectation maximization algorithm}

\title{
Multiple Data Mixture Models
}

\description{
Fits a layered or chained mixture model to a list representing multiple sources of data, using a choice of distributions and number of components for each data source.
}

\usage{
mdmixmod(X, K, K0=min(K), topology=LC_TOPOLOGY, 
         family=NULL, iter.max=LC_ITER_MAX)
\method{print}{mdmixmod}(mod, ...)
}

\arguments{
\item{X}{a list of observed data sources; the elements must be numeric vectors, matrices, or data frames.  Each element of \code{X} must have the same data length, that is, length of vectors or number of rows of matrices and data frames; but each element of \code{X} may be of arbitrary width, that is, number of columns of matrices and data frames.}
\item{K}{the vector of numbers of mixture components for the hidden variables corresponding to each observed data source.  If \code{length(K) < length(X)}, the values in \code{K} will be recycled; thus \code{K} may be a scalar, in which case the same number of components will be used throughout.  If \code{length(K) > length(X)}, then \code{K} will be truncated and a warning will be given.}
\item{K0}{the number of mixture components for the top-level hidden variable.}
\item{topology}{one of the model topologies in \code{\link{LC_TOPOLOGY}}, either \code{"layered"} or \code{"chained"}.  By default, \code{"layered"} is used.  Partial matches are allowed.}
\item{family}{a vector of names of distribution families to be used in fitting the models for each observed data source; each element of \code{family} must be one of the distribution names in \code{\link{LC_FAMILY}}.  The usual recycling and truncation rules are followed.  By default, \code{"normal"} is used and recycled to \code{c("normal", \dots, "normal")}.  Partial matches are allowed.}
\item{iter.max}{the maximum number of iterations for the EM algorithm, by default equal to \code{\link{LC_ITER_MAX}}.}
\item{mod}{an object of class \code{mdmixmod}.}
\item{...}{further arguments to \code{print}.}
}

\details{
In the layered model, a top-level hidden categorical random variable \eqn{\mathcal{Y}_0}{Y_0}, which can take on values from 1 to some positive integer \eqn{K_0}, generates categorical hidden random variables \eqn{\mathcal{Y}_1, \dots, \mathcal{Y}_Z}{Y_1, \dots, Y_Z} for some positive integer \eqn{Z}.  For \eqn{z = 1,\dots,Z}, each \eqn{\mathcal{Y}_z}{Y_z} can take on values from 1 to some positive integer \eqn{K_z}.  In the chained model, \eqn{\mathcal{Y}_0}{Y_0} generates \eqn{\mathcal{Y}_1}{Y_1}, which in turn generates \eqn{\mathcal{Y}_2}{Y_2}, etc., up to \eqn{\mathcal{Y}_{Z-1}}{Y_{Z-1}}, which generates \eqn{\mathcal{Y}_Z}{Y_Z}.

In both models, the \eqn{\mathcal{Y}_z}{Y_z}'s generate the observed mixture random variables \eqn{\mathcal{X}_1, \dots, \mathcal{X}_Z}{X_1, \dots, X_Z}, from which the elements of the observed data \code{X} are assumed to be drawn.  (That is, \code{Z = length(X)}, the number of list elements in \code{X}.)  The relationship between each \eqn{\mathcal{Y}_z}{Y_z} and \eqn{\mathcal{X}_z}{X_z} is the same as the relationship between \eqn{\mathcal{Y}}{Y} and \eqn{\mathcal{X}}{X} in \code{\link{mixmod}}.

As in \code{\link{mixmod}}, the EM algorithm attempts to maximize the Q-value, that is, the expected complete data (hidden and observed variables) log-likelihood.
}

\value{
A list of class \code{mdmixmod}, a subclass of \code{mixmod}, having the following elements:
\item{N}{the length of the data, that is, \code{length(X[[1]])} if \code{X[[1]]} is a vector, or \code{nrow(X[[1]])} if \code{X[[1]]} is a matrix or data frame.}
\item{Z}{the size of the data, that is, \code{Z = length(X)}.}
\item{D}{the vector of widths of the data, that is, \code{D[z] = 1} if \code{X[[z]]} is a vector, or \code{D[z] = ncol(X[[z]])} if \code{X[[z]]} is a matrix or data frame.}
\item{K}{the vector of the numbers of components in the lower-level mixture models.}
\item{K0}{the number of components in the top-level mixture model, that is, \eqn{K_0}{K_0}.}
\item{X}{the original data, with data frames converted to matrices.  If the elements of \code{X} were not named, they will be named \code{"X1",\dots,"XZ"} here.}
\item{npar}{the total number of parameters in the model.}
\item{npar.hidden}{the number of parameters for the hidden component portion of the model.}
\item{npar.observed}{the number of parameters for the observed data portion of the model.}
\item{iter}{the number of iterations required to fit the model.}
\item{params}{the parameters estimated for the model.  This is a list with elements \code{hidden} and \code{observed}, corresponding to distribution for the hidden and observed portions of the model.  \code{hidden} has elements \code{prob0}, the vector of probabilities for the possible values for \eqn{\mathcal{Y}_0}{Y_0}, and \code{cprob}, the list of matrices of conditional probabilities for the possible values of the \eqn{\mathcal{Y}_z}{Y_z}'s.  In the layered model, these are \eqn{K_0 \times K_z}{K_0-by-K_z} matrices of which the \eqn{(k_0,k_z)}th element represents \eqn{P(\mathcal{Y}_z = k_z \mid \mathcal{Y}_0 = k_0)}{P(Y_z = k_z | Y_0 = k_0)}.  In the chained model, the \eqn{z}th element of \code{cprob} is a \eqn{K_{z-1} \times K_z}{K_{z-1}-by-K_z} matrix of which the \eqn{(k_{z-1},k_z)}th element represents \eqn{P(\mathcal{Y}_z = k_z \mid \mathcal{Y}_{z-1} = k_{z-1})}{P(Y_z = k_z | Y_{z-1} = k_{z-1})}.  The chained model \code{hidden} also has elements \code{probz}, a list of vectors representing the marginal probabilities for the \eqn{\mathcal{Y}_z}{Y_z}'s, and \code{rprob}, representing the conditional probabilities of the \eqn{\mathcal{Y}_z}{Y_z}'s given the values of the \eqn{\mathcal{Y}_{z+1}}{Y_{z+1}}'s.  The elements of \code{observed} depend on the distribution family chosen in fitting the model.}
\item{stats}{a vector with named elements corresponding to the number of iterations, log-likelihood, Q-value, and BIC for the estimated parameters.}
\item{weights}{a list with elements \code{U}, a \eqn{N \times K_0}{N-by-K_0} matrix; \code{V}, a \eqn{Z}-length list of \eqn{K_0}-length lists of \eqn{N \times K_z}{N-by-K_z} matrices in the layered model, or a \eqn{Z}-length list of \eqn{K_{z-1}}-length lists of \eqn{N \times K_z}{N-by-K_z} matrices in the layered model; and \code{W}, a \eqn{Z}-length list of \eqn{N \times K_z}{N-by-K_z} matrices; representing weights used in the M-step of the EM algorithm for estimating the final set of parameters for the observed data portion of the model.}
\item{pdfs}{a list with elements \code{G}, \code{alpha}, \code{beta}, and \code{gamma}, representing various estimated density functions for the data.  \code{gamma} represents the estimated density of the observed data across all data sources under the fitted mixture model.}
\item{posterior}{the \eqn{N \times K_0}{N-by-K_0} matrix of which the \eqn{(n,k_0)}th element is the estimated posterior probability that the \eqn{n}th observation (across all data sources) was generated by the \eqn{k_0}th component.  Equal to the \code{U} element of \code{weights}.}
\item{assignment}{the vector of length \eqn{N} of which the \eqn{n}th element is the most probable top-level component to have generated the \eqn{n}th observation.  In other words, \code{assignment[n] = which.max(posterior[n,])}.}
\item{iteration.params}{a list of length \code{iter} giving the estimated parameters at each iteration of the algorithm.}
\item{iteration.stats}{a data frame of \code{iter} rows giving iteration statistics, as in \code{stats}, at each iteration of the algorithm.}
\item{family}{the vector of names of the distribution families used in the model. See \code{\link{LC_FAMILY}}.}
\item{distn}{the vector of names of the actual distributions used in the model. See \code{\link{LC_FAMILY}}.}
\item{iter.max}{the maximum number of distributions allowed in model fitting.}
\item{dname}{the name of the data.}
\item{dattr}{attributes of the data, used by model reporting functions to determine if the data have been scaled or otherwise transformed.}
}

\references{
McLachlan, G.J. and Thriyambakam, K.  (2008)  \emph{The EM Algorithm and Extensions}, John Wiley & Sons.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{LC_FAMILY}} for distributions and families; \code{\link{mixmod}} for fitting single-data mixture models; \code{\link{reporting_functions}} and \code{\link{qfun}} for reporting on objects of class \code{mdmixmod}; \code{\link{multiROC}} for performance evaluation; \code{convergence.plot} for behavior of the algorithm; \code{\link{simulation}} for simulating from the parameters of a model.
}

\examples{\dontrun{
data(exampleData3)
mod3 <- mdmixmod(exampleData3$X, c(2,3), 2)
    # univariate/multivariate normal, layered model
mod3$stats
#   iter       llik       qval        bic 
# 35.000  -9194.691 -12900.850 -18755.493

data(exampleData4)
mod4 <- mdmixmod(exampleData4$X, c(2,3), 2, family=c("wei", "gam"), 
                 topology="chained")
    # univariate Weibull (weisd) / multivariate gamma, chained model
mod4$stats
#    iter      llik      qval       bic 
# 121.000 -3451.138 -7527.333 -7268.387
}}

\keyword{graphs}
\keyword{models}
\keyword{cluster}
\keyword{classif}
