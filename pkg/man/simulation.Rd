\name{simulation}
\alias{simulation}
\alias{simulate.mixdata}
\alias{simulate.mdmixdata}
\alias{simulate.from.fit}
\alias{simulate.from.fit.mixmod}
\alias{simulate.from.fit.mdmixmod}

\title{
Generate Simulated Data
}

\description{
Simulation functions allow the generation of simulated data suitable for model fitting with \code{\link{mixmod}} or \code{\link{mdmixmod}} from either a \emph{de novo} set of parameters or the parameters of a fitted model.
}

\usage{
simulate.mixdata(N, distn, params)
simulate.mdmixdata(N, distn, params, topology=LC_TOPOLOGY)
simulate.from.fit(mod)
}

\arguments{
\item{N}{the length of the data to be simulated, that is, the number of independent simulated observations.}
\item{distn}{the name of the distribution for \code{simulate.mixdata}, or the vector of names of distributions for \code{simulate.mdmixdata}, to be used in generating the observed portion of the data.  Must be one of the distribution names (\emph{not} family names) in \code{\link{LC_FAMILY}}.}
\item{params}{the parameters of the complete data (hidden and observed) distributions; a list formatted as the \code{params} element of an object of class \code{mixmod} or \code{mdmixmod}.}
\item{topology}{one of the model topologies in \code{\link{LC_TOPOLOGY}},  either \code{"layered"} or \code{"chained"}.}
\item{mod}{an object of class \code{mixmod} or \code{mdmixmod}.}
}

\details{
\code{simulate.mixdata} and \code{simulate.mdmixdata} are for generating data from \emph{de novo} parameters, while \code{simulate.from.fit} generates data from the parameters of a fitted model.  The \code{params} argument must be a list formatted in the same way as the \code{params} element of a \code{mixmod} or \code{mdmixmod} object, such as is returned by the function \code{\link{mixmod}} or \code{\link{mdmixmod}}, respectively.  However, in the case of simulating a multiple data model from a layered topology, that is, \code{simulate.mdmixdata(..., topology="layered")}, the \code{probz} and \code{rprob} elements of \code{params$hidden} are not required.
}

\value{
See the Format section of \code{\link{exampleData}}.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}} for specifications of the parameters; \code{\link{exampleData}} for examples of the format of the data returned by simulation functions; \code{\link{constants}} for distribution names and topologies; \code{\link{mvnorm}}, \code{\link{mvweisd}}, \code{\link{weisd}}, \code{\link{mvgamma}} for relevant distributions; \code{\link{multiROC}} for performance evaluations using complete data.
}

\examples{\dontrun{
### create a layered univariate/multivariate normal data set,
### "edata", similar to exampleData3, from de novo parameters

# ensure reproducibility

set.seed(123)

# top-level (Y_0) probabilities; K_0 = 2

p0 = c(0.3, 0.7)

# second-layer conditional probabilities; K_1 = 2, K_2 = 3

Q1 = matrix(c(0.9, 0.1,
              0.3, 0.7),
            nrow=2, byrow=TRUE)

Q2 = matrix(c(0.6, 0.2, 0.2,
              0.1, 0.8, 0.1),
            nrow=2, byrow=TRUE)

Q = list(edata1=Q1, edata2=Q2)

# observed data parameters conditional on Y_1, Y_2

mu1 = c(3, -3)
mu2 = matrix(c( 3,  2,  1,  4,
                1, -2,  0, -1,
               -2, -4, -3, -5),
             nrow=3, byrow=TRUE)
colnames(mu2) = paste("edata2", 1:4, sep="")

sigmasq1 = 1 + rexp(2)
sigma1 = sqrt(sigmasq1)

Sigma2 = lapply(1:3, function(k) {
    Z = matrix(rnorm(100), ncol=4)
    colnames(Z) = paste("edata2", 1:4, sep="")
    mlecov(Z)
})

# put all parameters together into appropriate format

theta1 = list(mean=mu1, var=sigmasq1, sd=sigma1)

theta2 = list(mean=mu2, cov=Sigma2)

eparams = list(hidden   = list(prob0=p0, cprob=Q),
               observed = list(edata1=theta1, edata2=theta2))

# simulate

edata = simulate.mdmixdata(1000, c("norm", "mvnorm"), eparams)

# fit a model

edata.fit = mdmixmod(edata$X, c(2,3))

edata.fit
# Layered (normal, normal) mixture model ('norm', 'mvnorm')
# Data 'edata$X' of size 1000-by-(1,4) fitted to 2 (2,3) components
# Model statistics:
#      iter      llik      qval       bic 
#     10.00  -8420.70 -12108.18 -17207.51


### simulate from the fitted model and fit again

edata.fit.sim = simulate.from.fit(edata.fit)

edata.fit.sim.fit = mdmixmod(edata.fit.sim$X, c(2,3))

edata.fit.sim.fit
# Layered (normal, normal) mixture model ('norm', 'mvnorm')
# Data 'edata.fit.sim$X' of size 1000-by-(1,4) fitted to 2 (2,3) components
# Model statistics:
#       iter       llik       qval        bic 
#     11.000  -8419.721 -12124.658 -17205.554 
    }
}

\keyword{datasets}
\keyword{datagen}
