\name{multiROC}

\alias{multiROC}
\alias{plot.multiROC}
\alias{plot.multiROC.bw}

\concept{ROC curve}

\title{
ROC Curves for Multiple Predictions
}

\description{
Functions to calculate and plot ROC (Receiver Operating Characteristic) curves for one or more sets of predictions.  All functions listed here require the suggested package \code{ROCR}.
}

\usage{
multiROC(predictions, labels, quasi=FALSE)

\method{plot}{multiROC}(mroc, do.legend=TRUE, do.auc=TRUE, do.dca=TRUE,
    xlab=ifelse(mroc$quasi, "all positive rate", "false positive rate"),
    ylab="true positive rate", main="", colorplus=0, 
    colors=colorplus+(1:length(mroc$perf)), ...)

plot.multiROC.bw(mroc, do.legend=TRUE, do.auc=TRUE, do.dca=TRUE,
    xlab=ifelse(mroc$quasi, "all positive rate", "false positive rate"),
    ylab="true positive rate", main="", ltys=1:length(mroc$perf), ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{predictions}{A list with named elements, each element being a vector of predictions, that is, numeric values between 0 and 1 inclusive.  All vectors must be of the same length, say \eqn{N}.  May be a data frame.}
\item{labels}{a vector of logical values of length \eqn{N}, in which a value of \code{TRUE} indicates a true positive; or a list of such vectors of the same length, with the same element names, as \code{predictions}.  May be a data frame.}
\item{quasi}{logical; if \code{TRUE}, the x-axis of the ROC curve denotes the total positive rate rather than the false positive rate.}
\item{mroc}{an object of class \code{multiROC}.}
\item{do.legend}{logical; if \code{TRUE}, a legend for the various ROC curves will be shown on the plot.}
\item{do.auc}{logical; if \code{TRUE}, AUC, or Area Under the Curve, will be shown in the legend.  Larger AUC indicates better performance.  No effect if \code{do.legend} is \code{FALSE}.}
\item{do.dca}{logical; if \code{TRUE}, DCA, or Distance of Closest Approach of the curve to (0,1), will be shown in the legend.  Smaller DCA indicates better performance.  No effect if \code{do.legend} is \code{FALSE}.}
\item{xlab}{title for the x-axis.}
\item{ylab}{title for the y-axis.}
\item{main}{overall plot title.}
\item{colorplus}{offset for colors used in plotting curves; colors start at \code{colors[[1]]+colorplus}.}
\item{colors}{colors used in plotting curves, by default \code{1, 2, \dots, length(predictions)}, as in the code{col} argument to \code{plot} in package \code{graphics}.}
\item{ltys}{line types used in plotting curves.}
\item{...}{further arguments to \code{plot}.}
}

\details{
\code{multiROC} is used to evaluate the performance of one or more sets of predictors against either a single set of true-positive / true-negative labels, or against a list of such sets of the same length as the list of predictions.  A quasi-ROC curve, that is, \code{multiROC(\dots, quasi=TRUE)}, is used when some number of true positives is known, but true negatives are not; that is, it is reasonable to assume that there exists some number of true positives among the \code{FALSE} labels, so long as the number of these ``unknown positives'' is small relative to the total number of samples.

See package \code{ROCR} for more information on ROC curves and other methods of predictor evaluation.
}

\value{
The \code{plot} methods are used for their side-effects; \code{plot.multiROC.bw} plots in black and white.  For \code{multiROC}, an object of class \code{multiROC}, a list with the following elements:
\item{pred}{a list of objects of class \code{prediction}.  See \code{prediction-class} in package \code{ROCR}.}
\item{perf}{a list of objects of class \code{performance}.  See \code{performance-class} in package \code{ROCR}.}
\item{meas}{a list of vectors, each vector having named elements \code{auc} and \code{dca}, giving the AUC and DCA for each prediction vector in \code{predictions}.}
\item{quasi}{logical; same as the \code{quasi} argument in the function call.}
}

\author{
Daniel Dvorkin
}

\note{
In future releases, \code{plot.multiROC.bw} may be deprecated, with black-and-white plotting being made an option in the main \code{plot} method.
}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}}, \code{\link{exampleData}}, \code{\link{simulation}}, \code{\link{convergence.plot}}; package \code{ROCR}; \code{plot} in package \code{graphics}.
}

\examples{\dontrun{
## multiple predictors, single label vector:
## demonstrate that combining data in the layered model is
## more effective than using individual data sources

data(exampleData3)
mod.ed31 = mixmod(exampleData3$X$ed31, 2)
mod.ed32 = mixmod(exampleData3$X$ed32, 3)
mod.layered = mdmixmod(exampleData3$X, c(2,3), 2)
plot(multiROC(list(ed31=mod.ed31$posterior[,1],
                   ed32=mod.ed32$posterior[,1],
                   layered=mod.layered$posterior[,1]),
              (exampleData3$Y0 == 1)))

## multiple predictors, multiple labels:
## compare univariate normal 2-component normal 
## with multivariate 3-component Weibull

data(exampleData1)
mod1 = mixmod(exampleData1$X, 2)
data(exampleData2)
mod2 = mixmod(exampleData2$X, 3, family="wei")
plot(multiROC(list(mod1=mod1$posterior[,1],
                   mod2=mod2$posterior[,1]),
              list(mod1=(exampleData1$Y==1),
                   mod2=(exampleData2$Y==1))))
}}

\keyword{hplot}
