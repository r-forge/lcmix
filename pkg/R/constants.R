### constants.R:  "constants" for the lcmix package

LC_LOG2 <- log(2)
LC_LOGPI <- log(pi)
LC_LOG2PI <- log(2*pi)
LC_LOGSQRT2PI <- log(sqrt(2*pi))
LC_LOGSQRT2 <- log(sqrt(2))
	# it's just useful to have these around

LC_EPS   <- sqrt(.Machine$double.eps)
LC_LEPS  <- log(LC_EPS)
LC_REPS  <- 1.0 - LC_EPS
LC_LREPS <- log(LC_REPS)
	# LC_EPS and LC_REPS ("reverse LC_EPS") are numbers very close, but not quite equal, to 0 and 1 respectively

LC_MAXSTDNORM <- 8
LC_MINSTDNORM <- -LC_MAXSTDNORM
	# maximum and minimum values of standard normal random variables; going outside these bounds with a variable that is supposed to be standard normal (e.g. Phi^-1(F(x)) for any random variable x) may lead to instability

LC_ITER_MAX <- 1000
	# default max iterations for EM
LC_ITER_TOL <- 1e-6 # 1e-6 or LC_EPS
	# point below which relative changes in log-likelihood must fall to terminate the algorithm

LC_TOPOLOGY <- c("layered", "chained")
	# topologies for multi-data models

LC_FAMILY <- list(normal      = list(uni="norm",  multi="mvnorm"),
                  pvii        = list(uni="pvii",  multi="mvpvii"),
                  weibull     = list(uni="weisd", multi="mvweisd"),
	              gamma       = list(uni="gamma", multi="mvgamma"),
	              exponential = list(uni="exp",   multi="mvexp"),
	              altpvii     = list(uni="apvii", multi="amvpvii"))
	# names of acceptable model distribution _families_ are names(LC_FAMILY); the $uni and $multi elements of the elements of LC_FAMILY give the names of the applicable univariate and multivariate distributions, respectively, for those families

LC_NONNEGFAM <- c("weibull", "gamma", "exponential")
	# distribution families requiring non-negative values

LC_DISTN <- as.vector(unlist(LC_FAMILY))
	# names of acceptable model _distributions_


LC_SIMPAR = list(
	
	hidden = list(
	
		prob0 = c(0.03, 0.97),
		
		probz = list(
			binding      = c(0.03, 0.97),
			expression   = c(0.05, 0.92, 0.03),
			conservation = c(0.25, 0.75)
		)
	),
	
	observed = list(
	
		binding = list(
			mean = c(0.825, -0.825),
			var  = c(1.50, 0.50),
			sd   = sqrt(c(1.50, 0.50))
		),
		
		expression = list(
			
			mean = matrix(c( 1.5,  1.5,  1.5,
			                 0.00,  0.00,  0.00,
			                -1.5, -1.5, -1.5),
			              nrow=3, byrow=TRUE),
			
			cov = list(
				matrix(c( 3.00,  1.00, -0.50,
				          1.00,  3.00,  0.00,
				         -0.50,  0.00,  3.00),
				        nrow=3, byrow=TRUE),
				diag(3),
				matrix(c( 2.00, -0.20,  0.10,
				         -0.20,  2.00,  0.10,
				          0.10,  0.10,  2.00),
				        nrow=3, byrow=TRUE)
			)
		),
		
		conservation = list(
			mean = c(1.25, -1.25),
			var  = c(1.50, 0.50),
			sd   = sqrt(c(1.50, 0.50))
		)
	)	
) # marginal distribution parameters for simulating data with similar performance to the Ci data