### constants.R:  "constants" for the lcmix package

LC_LOG2PI <- log(2*pi)
	# it's just useful to have this around

LC_EPS   <- sqrt(.Machine$double.eps)
LC_LEPS  <- log(LC_EPS)
LC_REPS  <- 1.0 - LC_EPS
LC_LREPS <- log(LC_REPS)
	# LC_EPS and LC_REPS ("reverse LC_EPS") are numbers very close, but not quite equal, to 0 and 1 respectively

LC_MAXSTDNORM <- 8
LC_MINSTDNORM <- -LC_MAXSTDNORM
	# maximum and minimum values of standard normal random variables; going outside these bounds with a variable that is supposed to be standard normal (e.g. Phi^-1(F(x)) for any random variable x) may lead to instability

LC_EPSEPS <- 1e-300
	# for when LC_EPS just isn't small enough; use with caution!

LC_ITER_MAX <- 10000
	# default max iterations for EM

LC_TOPOLOGY <- c("layered", "chained")
	# topologies for multi-data models

LC_FAMILY <- list(normal  = list(uni="norm",  multi="mvnorm"),
                  weibull = list(uni="weisd", multi="mvweisd"),
	              gamma   = list(uni="gamma", multi="mvgamma"))
	# names of acceptable model distribution families are names(LC_FAMILY); the $uni and $multi elements of the elements of LC_FAMILY give the names of the applicable univariate and multivariate distributions, respectively, for those families

LC_NONNEGFAM <- c("weibull", "gamma")
	# distribution families requiring non-negative values
