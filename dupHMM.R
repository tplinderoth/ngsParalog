#!/usr/bin/env Rscript

# dupHMM.R

# TODO:
###############################################################################

library(expm, quietly=TRUE, warn.conflicts=FALSE)
library(truncnorm, quietly=TRUE)
library(docopt, quietly=TRUE)

###### functions ######

initializePi <- function() {
	# initialize initial distribution vector uniformly
	
	pi <- c(0.5, 0.5)
	return(pi)
}

initializeP <- function() {
	# initialize transition probability matrix uniformly
	
	nstates <- 2
	trans <- matrix(NA, nrow=nstates, ncol=nstates)

	p <- 1/nstates
	for (i in 1:nstates) {
		for (j in 1:nstates) {
			trans[i,j] <- p
		}
	}
	
	return(trans)
}

fitlr <- function (lr, tailcutoff = 1.0) {
	# lr: likelihood ratios
	# tailcutoff: fit alternative distribution to LRs under the 'tailcutoff' quantile (avoid influence of extreme values)
	
	# check params
	if (tailcutoff <= 0 || tailcutoff > 1) stop("LR quantile cutoff", tailcutoff, "outside of (0,1]")
	
	# truncate sites with extremeley high LR
	truncidx <- which(lr < quantile(lr, tailcutoff))
	sublr <- lr[truncidx]
	
	# estimate proportion of LRs under the null
	lrquantile <- 5.411894 # 99% quantile for null LR distribution
	dupidx <- which(sublr > lrquantile)
	if (length(dupidx) == 0) dupidx <- which(sublr > 0)
	
	pnull <- 1.0 - length(dupidx)/length(sublr)
	pmax <- 0.999999
	pmin <- 1-pmax
	if (pnull > pmax) pnull <- pmax
	
	# approximate mean of alternative LR values
	ncpguess <- mean(sublr[dupidx])
	ncpmax <- max(sublr)
	
	fit <- tryCatch(
			optim(par=c(pnull, ncpguess), fn=lrllh, method="L-BFGS-B", lower=c(pmin, 0), upper=c(pmax, ncpmax), x=sublr)$par,
			error = function(err) {
				cat(paste("LR DISTRIBUTION OPTIMIZATION FAILURE:\n", err, sep=''))
				return(NULL)
			})
	
	return(fit)
}

lrllh <- function (par, x) {
	# LR ~ par[1]*[0.5*chisq(df=0, ncp=0) + 0.5*chisq(df=1, ncp=0)] + (1-par[1])*chisq(df=1, ncp=par[2])
	# par[1] = probability LR comes from null distribution 
	# par[2] = noncentrality parameter for alternate chisquare distribution
	
#	llh <- rep(NA, length(x))
#	for (i in 1:length(x)) {
#		if (x[i] == 0) llh[i] <- par[1]*0.5 else llh[i] <- par[1]*0.5*dchisq(x[i], df=1, ncp=0) + (1-par[1])*dchisq(x[i], df=par[2], ncp=par[3])
#	}
#	llh <- replace(llh, llh==0, .Machine$double.xmin)
#	ret <- -sum(log(llh))
	
	# faster in R:
	zeroidx <- which(x==0)
	nzero <- length(zeroidx)
	if (nzero != 0) x <- x[-zeroidx]
	
	dzero <- par[1]*0.5
	y1 <- rep(dzero, nzero)
	y2 <- dzero*dchisq(x, df=1, ncp=0) + (1-par[1])*dchisq(x, df=1, ncp=par[2])
	
	llh <- c(y1, y2)
	llh <- replace(llh, llh==0, .Machine$double.xmin)
	-sum(log(llh))
}

covllh <- function (par, x) {
	# coverage ~ par[1]*truncnormal(mean=par[2], sd=par[3], a=0, b=Inf) + (1-par[1])*truncnorm(mean=2*par[2], sd=par[4], a=par[5], b=Inf)
	# par[1] = probability x is from null distribution
	# par[2] = null truncated normal mean
	# par[3] = null truncated normal standard deviation
	# par[4] = alternative truncated normal standard deviation
	# par[5] = alternative truncated normal lower bound
	
	altmean <- 2*par[2]
	if (par[5] > altmean) par[5] <- altmean # to avoid NaN
	llh <- par[1] * dtruncnorm(x, mean=par[2], sd=par[3], a=0, b=Inf) + (1-par[1]) * dtruncnorm(x, mean=altmean, sd=par[4], a=par[5], b=Inf)
	llh <- replace(llh, llh==0, .Machine$double.xmin)
	-sum(log(llh))
}

fitCoverage <- function(coverage, lr, maxcov=NULL, alt_min=0) {
	# coverage: vector of average individual coverage for each site
	# lr: duplication likelihood ratios
	# maxcov: truncate coverage above maxcov when fitting distribution
	# alt_min: truncated normal duplicated coverage distribution lower bound
	
	# truncate coverage
	if (!is.null(maxcov)) {
		if (maxcov <= 0) stop("Maximum coverage has to be > zero")
		if (alt_min > maxcov) stop(paste("--dupcovmin", alt_min, "greater than --maxcoverage", maxcov, "\n"))
		keepidx <- which(coverage <= maxcov)
		coverage <- coverage[keepidx]
		lr <- lr[keepidx]
	}
	
	# approximate mean and variance of average individual coverage
	nullcutoff <- min(lr)
	nondupidx <- which(lr <= nullcutoff)
	inc <- 0.5
	while(length(nondupidx) < 2) {
		nullcutoff <- nullcutoff + inc
		nondupidx <- which(lr <= nullcutoff)
	}
	
	nullavg <- mean(coverage[nondupidx])
	nullsd <- sqrt(var(coverage[nondupidx]))
	null_mu_min <- min(coverage)
	null_mu_max <- 3*nullavg
	if (!is.null(maxcov) && null_mu_max > maxcov) null_mean_max <- 0.5 * max(coverage)
	
	alt_max <- 2*null_mu_max # maximum value for alternative truncated normal lower bound	
	if (alt_min > alt_max) {
		cat(paste("--dupcovmin", alt_min, "is greater than", alt_max, "(maximum value for duplicated distribution lower bound)\n"))
		return(NULL)
	}
	
	# approximate proportion of nonduplicated sites
	pnullcutoff <- 5.411894 # 99% quantile for null LR distribution
	pnull <- length(which(lr < pnullcutoff))/length(lr)
	pmax <- 0.999999
	pmin <- 1 - pmax
	if (pnull > pmax) pnull <- pmax

	# set starting point for optimization
	
	altsd <- 2*nullsd
	alt_lowbound <- alt_min
	
	if (!is.null(maxcov)) {
		if (alt_max > maxcov) cat(paste("WARNING: Upper bound for duplicated mean coverage of", alt_max, "exceeds --maxcoverage", maxcov, "\n"))
		if (alt_lowbound > maxcov) cat(paste("WARNING: Duplicated coverage lower bound guess of", alt_lowbound, "exceeds --maxcoverage", maxcov, "\n"))
	}
	
	# find MLE parameter values
	fit <- tryCatch( 
			optim(par=c(pnull, nullavg, nullsd, altsd, alt_lowbound), fn=covllh, method="L-BFGS-B", 
					lower=c(pmin, null_mu_min, 1e-12, 1e-12, alt_min), upper=c(pmax, null_mu_max, Inf, Inf, alt_max), x=coverage)$par,
			error = function(err) {
				print(paste("COVERAGE DISTRIBUTION OPTIMIZATION FAILURE:", err, sep=''))
				return(NULL)
			})
	
	return(fit)
}

penalizeLR <- function(lr, penalty, sampn=NULL) {
	# lr: vector of duplication likelihood ratios
	# penalty: type of penalized likelihood ratios to use
	# sampn: sample size (number of diploid individuals)
	
	if (penalty %in% c('aic', 'AIC')) {
		lr <- lr - 2
	} else if (penalty %in% c('bic', 'BIC')) {
		if (is.null(sampn)) stop("Use of BIC penalized likelihood ratio requires sample size")
		adjust <- log(sampn)
		lr <- lr - adjust
	} else stop(paste(penalty, "is an invalid likelihood penalization"))
	lr[lr < 0] <- 0
	
	return(lr)
}

fitEmissions <- function(lr, coverage=NULL, emittype, lrmax_quantile, maxcoverage=NULL, min_alt_cov=0) {
	# estimates parameters used to calculate emission probabilities
	
	# lr: vector of duplicaiton likelihood ratios
	# coverage: vector of average individual coverage
	# emittype: use (0) LR only or (1) LR and coverage as emissions for hmm
	# lrmax_quantile: ignore LRs above this quantile
	# maxcoverage: maximum coverage
	# min_alt_coverage: lower bound for truncated normal minimum parameter of alternative coverage distribution

	# check function params
	if (! emittype %in% c(0,1)) stop(paste("Invalid emission type", emittype))
	
	# estimate parameters for null coverage distribution
	covparams <- NULL
	if (emittype == 1) {
		cat("\nFitting coverage distribution\n")
		covpar <- fitCoverage(coverage=coverage, lr=lr, maxcov=maxcoverage, alt_min=min_alt_cov)
		if (is.null(covpar)) stop("Unable to fit coverage distribution")
		
		# print fitted coverage distribution params
		covparams <- c(covpar[1], covpar[2], covpar[3], 2*covpar[2], covpar[4], covpar[5]) # assume duplicated sites have 2x the coverage as normal sites
		names(covparams) <- c("cov_ndProb", "cov_ndMean", "cov_ndSD", "cov_dupMean", "cov_dupSD", "cov_dupMin")
		print(covparams)
	}
	
	# estimate parameters for LR distribution
	cat("\nFitting LR distribution\n")
	lrpar <- fitlr(lr=lr, tailcutoff=lrmax_quantile)
	if (is.null(lrpar)) stop("Unable to fit LR distribution")
	
	# print fitted LR distribution params
	lrparams <- c(lrpar[1], lrpar[2])
	names(lrparams) <- c("lr_ndProb", "lr_dupNCP")
	print(lrparams)
	
	# return parameter estimates
	fits <- list(lrparams, covparams)
	
	return(fits)
}

initializeEmissions <- function(lr, lrpar, coverage=NULL, covpar=NULL, emittype) {
	# initialize emission probabilites
	
	# lr: vector of duplicaiton likelihood ratios
	# larpar: vector with LR density params: [1] nondup prob, [2] dup noncentrality param
	# coverage: vector of average individual coverage
	# covpar: vecotr with coverage density params: [1] nondup prob, [2] nondup mean, [3] nondup std deviation, [4] dup mean, [5] dup std deviation, [6] dup min
	# emittype: use (0) LR only or (1) LR and coverage as emissions for hmm
	
	# check function params
	if (! emittype %in% c(0,1)) stop(paste("Invalid emission type", emittype))

	# calculate emission probabilities
	cat("\nCalculating emission probabilities\n")
	
	# calculate probabilites of the observed LR
	lrseq <- seq(from=2, to=max(ceiling(lr)), by=1)
	lrprob <- matrix(NA, nrow=2, ncol=length(lrseq)+1)
	
	lrmax <- 1
	seenalt <- 0
	lrprob[1,1] <- 0.5 + 0.5*pchisq(1,1)
	lrprob[2,1] <- pchisq(1, df=1, ncp=lrpar[2])
	nullsum <- lrprob[1,1]
	altsum <- lrprob[2,1]
	
	for (i in lrseq) {
		lrprob[1,i] <- (0.5 + 0.5*pchisq(q=i, df=1)) - nullsum # null
		palt <- pchisq(i, df=1, ncp=lrpar[2])
		if (palt >= altsum) lrprob[2,i] <- palt - altsum else break # alternative (need this to check for R numerical instability)
		nullsum <- nullsum + lrprob[1,i]
		altsum <- altsum + lrprob[2,i]
		if (seenalt == 0 && lrprob[2,i] > 0) seenalt <- 1 # check if the alternative distribution has been entered 
		if (lrprob[1,i] == 0 && lrprob[2,i] == 0 && seenalt == 1) break
		lrmax <- lrmax + 1
	}
	
	# calculate probabilities of the observed coverage
	covprob <- NULL

	covmax <- 1
	if (emittype == 1) {
		covseq <- seq(from=2, to=max(ceiling(coverage)), by=1)
		covprob <- matrix(NA, nrow=2, ncol=length(covseq)+1)
		
		seenalt <- 0
		covprob[1,1] <- ptruncnorm(1, mean=covpar[2], sd=covpar[3], a=0, b=Inf) # null
		covprob[2,1] <- ptruncnorm(1, mean=covpar[4], sd=covpar[5], a=covpar[6], b=Inf) #alternative
		nullsum <- covprob[1,1]
		altsum <- covprob[2,1]
	
		for (i in covseq) {
			covprob[1,i] <- ptruncnorm(i, mean=covpar[2], sd=covpar[3], a=0, b=Inf) - nullsum # null
			covprob[2,i] <- ptruncnorm(i, mean=covpar[4], sd=covpar[5], a=covpar[6], b=Inf) - altsum # alternate
			nullsum <- nullsum + covprob[1,i]
			altsum <- altsum + covprob[2,i]
			if (seenalt == 0 && covprob[2,i] > 0) seenalt <- 1 # check if alternative distribution has been entered
			if (covprob[1,i] == 0 && covprob[2,i] == 0 && seenalt == 1) break
			covmax <- covmax + 1
		}
	}
	
	# bound the maximum LR and coverage values
	lr[lr==0] <- 1
	lr <- ceiling(lr)
	lr <- replace(lr, lr > lrmax, lrmax)
	
	if (emittype == 1) {
		coverage[coverage==0] <- 1
		coverage <- ceiling(coverage)
		coverage <- replace(coverage, coverage > covmax, covmax)
	} else coverage <- rep(1,length(lr))
	
	# check for excess zero emission probabilities
	if (sum(lrprob[1,], na.rm=TRUE) == 0) return(list("All nonduplicated LR emission probabilities are zero!", list(lr, coverage)))
	if (sum(lrprob[2,], na.rm=TRUE) == 0) return(list("All duplicated LR emission probabilities are zero!", list(lr, coverage)))
	if (!is.null(covprob) && sum(covprob[1,], na.rm=TRUE) == 0) return(list("All nonduplicated coverage emission probabilities are zero!", list(lr, coverage)))
	if (!is.null(covprob) && sum(covprob[2,], na.rm=TRUE) == 0) return(list("All duplicated coverage emission probabilities are zero! Try increasing --dupcovmin.", list(lr, coverage)))

	# calculate emission probabilities
	b <- list(matrix(NA, nrow=lrmax, ncol=covmax), matrix(NA, nrow=lrmax, ncol=covmax))
	
	if (emittype == 1) {
		# calculate joint probabilities of LR and coverage
		for (i in 1:lrmax) {
			for (j in 1:covmax) {
				b[[1]][i,j] <- lrprob[1,i] * covprob[1,j] # null
				b[[2]][i,j] <- lrprob[2,i] * covprob[2,j] # alternative
			}
		}
	} else {
		# only use LR probabilities
		b[[1]][,1] <- lrprob[1,1:lrmax] # null
		b[[2]][,1] <- lrprob[2,1:lrmax] # alternative
	}
	
	# assign uniform emission probability to cases with zero probability under both the null and alternative
	zeroidx <- which(b[[1]] == 0 & b[[2]] == 0) # indices that are zero (or < machine precision) under null and alternative
	if (length(zeroidx) > 0) {
		uniprob <- 1/(nrow(b[[1]])*ncol(b[[1]]))
		b[[1]][zeroidx] <- uniprob
		b[[2]][zeroidx] <- uniprob
	}
	
	# make the emission probs sum to one
	b[[1]] <- b[[1]]/sum(b[[1]]) # null emissions must sum to 1
	b[[2]] <- b[[2]]/sum(b[[2]]) # alternative emissions must sum to 1
	
	if (min(c(min(b[[1]][which(b[[1]]>0)]), min(b[[2]][which(b[[2]]>0)]))) < .Machine$double.xmin) {
		cat("WARNING: Smallest nonzero emission probability less than machine precision - consider decreasing LR quantile\n")
	}
	
	return(list(b, list(lr, coverage)))
}

seqPMatrix <- function(p, sites, logscale=0) {
	# calculates a list of n-step transition matrices based on the distance between snps

	# p: transition probability matrix
	# sites: vector steps of the HMM (SNP positions)
	
	P <- list(p) # P[[1]] is a dummy
	if (logscale) P[[1]] <- log(P[[1]])
	
	for (i in 2:length(sites)) {
		P[[i]] <- p %^% (sites[i] - sites[i-1])
		if (logscale) P[[i]] <- log(P[[i]])
	}
	
	return(P)
}

hmmForward <- function (obs, pi, p, b) {
	# forward algorithm
	# use scaling of probabilities for underflow protection (Rabiner 1989 pg. 16)
	
	# obs: vector of observations
	# pi: initial distribution vector
	# p: transition probability matrices
	# b: emission probability matrix
	
	nstates <- length(pi)
	T <- length(obs[[1]])
	alpha <- matrix(data=NA, nrow=T, ncol=nstates) # alpha(t,i) = P(O1, O2, O3,..., Ot, qt=Si), forward variables
	c <- rep(NA, T) # scaling variables
	
	# initialization step - might need to adjust this and start at alpha(0,i) as in Durbin
	c[1] <- 0
	for (i in 1:nstates) {
		alpha[1, i] <- pi[i] * b[[i]][ obs[[1]][1], obs[[2]][1] ] # alpha(1,i) = pi_i * b_i(O1)
		c[1] <- c[1] + alpha[1,i]
	}
	if (c[1] == 0) stop("Scaling variable is zero")
	
	# scale the alpha[1,i]
	c[1] <- 1/c[1]
	for (i in 1:nstates) alpha[1, i] <- c[1] * alpha[1, i]

	# induction step
	for (t in 2:T) {
		c[t] <- 0
		
		# alpha(t+1, i) = sum_i=1_to_N( alpha(t,j) * p_ji * b_i(O_t+1) )
		for (i in 1:nstates) {
			alpha[t, i] <- 0
			for (j in 1:nstates) alpha[t, i] <- alpha[t, i] + alpha[t-1, j] * p[[t]][j, i]
			alpha[t, i] <- alpha[t, i] * b[[i]][ obs[[1]][t], obs[[2]][t] ]
			c[t] <- c[t] + alpha[t, i]
		}
		if (c[t] == 0) stop("Scaling variable is zero")
		
		# scale the forward variables
		c[t] <- 1/c[t]
		for (i in 1:nstates) alpha[t, i] <- c[t] * alpha[t, i]
	}
	
	# return list of scaled forward variables and scaling variables
	return (list(alpha, c))
}


hmmBackward <- function(obs, p, b, c) {
	# backward algorithm
	# using scaling of of probabilities for underflow protection
	
	# obs: vector of observations
	# p: list of transition probability matrices
	# b: emission probability matrix
	# scaling variables from the forward algorithm
	
	nstates <- nrow(p[[1]])
	T <- length(obs[[1]])
	beta <- matrix(data=NA, nrow=T, ncol=nstates) # beta(t,i) = P(O_t+1, O_t+2,..., O_T | q_t = S_i), backward variables
	
	# initializaiton
	for (i in 1:nstates) {
		beta[T, i] <- c[T]
	}
	
	# recursion
	for (t in (T-1):1) {
		for (i in 1:nstates) {
			beta[t, i] <- 0
			for (j in 1:nstates) beta[t, i] <- beta[t, i] + p[[t+1]][i, j] * b[[j]][ obs[[1]][t+1], obs[[2]][t+1] ] * beta[t+1, j]
			beta[t, i] <- c[t] * beta[t, i]
		}
	}
	
	# return matrix of backwards variables
	return(beta)
}

hmmStateVars <- function (obs, alpha, beta, p, b) {
	# probability that p_ij is used at position t, P(x_t = i, x_t+1 = j | x, lambda)
	
	# obs: vector of observations
	# alpha: forward variables
	# beta: backward variables
	# p: transition probability matrices
	# b: emission probability matrix
	
	T <- nrow(alpha)
	nstates <- ncol(alpha)
	gamma <- lapply(1:(T-1), matrix, data=NA, nrow=nstates, ncol=nstates) # P(x_t = i, x_t+1 = j | x, lambda)
	psi <- matrix(0, nrow=T, ncol=nstates) # P(x_t = i | x, lambda)
	
	for (t in 1:(T-1)) {
		
		# calculate denominator, P(X)
		px <- 0
		for (i in 1:nstates) {
			for (j in 1:nstates) px <- px + ( alpha[t,i] * p[[t+1]][i, j] * b[[j]][ obs[[1]][t+1], obs[[2]][t+1] ] * beta[t+1, j] )
		}
		if (px == 0) stop("Probability of data is zero ")
		
		
		for (i in 1:nstates) {
			psi[t, i] <- 0
			for (j in 1:nstates) {
				# probability that p_ij is used at position t
				gamma[[t]][i,j] <- (alpha[t, i] * p[[t+1]][i, j] * b[[j]][ obs[[1]][t+1], obs[[2]][t+1] ] * beta[t+1, j]) / px
				# probability that x_t = i
				psi[t, i] <- psi[t, i] + gamma[[t]][i,j]
			}
		}
	}
	
	px <- 0
	for (i in 1:nstates)
		px <- px + alpha[T, i]
	
	if (px == 0) stop("Probability of data is zero")
	
	for (i in 1:nstates)
		psi[T, i] <- alpha[T, i]/px
	
	return(list(gamma, psi))
}

hmmEstParam <- function(gamma, psi, estimate=c(1,2,3), obs=NA) {
	# Baum-Welch parameter re-estimation step 
	# obs required if estimating emission matrix
	
	# gamma: probility that p_ij is used at step t
	# psi: P(x_t) = i at each step t
	# estimate: what parameters to restimate, 1=>pi, 2=>transition probabilities, 3=>emission probabilities
	# obs: vector of potential observations
	
	# check input parameters
	if (! all(estimate %in% c(1,2,3))) stop("HMM parameters to estimate not in (1) initial distribution, (2) transition matrix, (3) emission matrix")

	nstates <- ncol(psi)
	T <- nrow(psi)
	params <- list(NA, NA, NA) # [[1]] = pi, [[2]] = transition matrix, [[3]] = emission matrix
	
	# estimate initial distribution
	if (1 %in% estimate) {
		params[[1]] <- rep(NA, nstates)
		for (i in 1:nstates) params[[1]][i] <- psi[1, i]
	}
	
	# estimate transition probability matrix
	if (2 %in% estimate) {
		params[[2]] <- matrix(NA, nrow=nstates, ncol=nstates)
		
		for (i in 1:nstates) {
			for (j in 1:nstates) {
				expij <- 0 # expected number of transitions from state i to j
				exi <- 0 # expected number of transitions from state i - can store this for estimating B
				for (t in 1:(T-1)) {
					expij <- expij + gamma[[t]][i,j]
					exi <- exi + psi[t,i]
				}
				if (exi == 0) stop("Attempt to divide by zero")
				params[[2]][i,j] <- expij/exi
			}
		}
	}
	
	if (3 %in% estimate) {
		# the following code for re-estimating the emission probabilities is untested
		params[[3]] <- list(matrix(0, nrow=nrow(obs[[1]]), ncol=ncol(obs[[1]])), matrix(0, nrow=nrow(obs[[1]]), ncol=ncol(obs[[1]])))
		# estimate emission probability matrix
	#	for (i in 1:nstates) {
	#		for (j in 1:ncol(params[[3]])) {
	#			exobsj = 0 # expected number of times in state i observing O_j
	#			exi = 0 # expected number of times in state i
	#			for (t in 1:T) {
	#				if (obs[t] == j) exobsj = exobsj + psi[t, i]
	#				exi = exi + psi[t, i]
	#			}
	#			params[[3]][i, j] = exobsj/exi
	#		}
	#	}
	
		# alternative estimate for emission prob matrix (avoids loop over observation types) - only works when observation index = observation symbol
		for (i in 1:nstates) {
			exi <- 0 # expected number of times in state i
			for (t in 1:T) {
				params[[3]][[i]][ obs[[1]][t], obs[[2]][t] ] <- params[[3]][[i]][ obs[[1]][t], obs[[2]][t] ] + psi[t, i]
				exi <- exi + psi[t, i]
			}
			if (exi == 0) stop("Attempt to divide by zero")
			for (j in 1:nrow(params[[3]][[1]])) {
				for (k in 1:ncol(params[[3]][[1]])) params[[3]][[i]][j,k] <- params[[3]][[i]][j,k]/exi
			}
		}
	}

	return(params)
}

hmmProbObs <- function(c) {
	# returns log(P(X))
	
	# c: scaling variables from forward algorithm
	
	pobs <- 0
	for (j in 1:length(c)) pobs <- pobs + log(c[j])
	
	return(-pobs)
}

hmmBaumWelch <- function(obs, steps, pi, p, b, maxiter=100, pdifflimit=1e-4, est=c(1,2,3)) {
	# Baum Welch algorithm
	
	# obs: observations
	# pi: initital probability distribution
	# p: transition probability matrix
	# b: emission probability matrix
	# maxiter: maximum number of EM iterations
	# logdiff: minimum difference in the log likelihoods
	# parameters to estimate: 1=>inititial distribution, 2=>transition matrix, 3=>emission probs
	
	# check input parameters
	if (maxiter < 1) stop("Max number of EM algorithm iterations set to < 1")
	if (pdifflimit < 0) stop("EM algorithm difference in log probility set to < 0")
	
	lambda <- list(pi, p, b)
	oldlogPX <- 0
	logPX <- 0
	iter <- 0
	
	while (iter < maxiter) {
		seqp <- seqPMatrix(p=lambda[[2]], sites=steps, logscale=0)
		
		# calculate forward variables
		fvar <- hmmForward(obs=obs, pi=lambda[[1]], p=seqp, b=lambda[[3]])

		# calculate backward variables
		bvar <- hmmBackward(obs=obs, p=seqp, b=lambda[[3]], c=fvar[[2]])

		# calculate probability of i->j and of being in state i at time t
		svar <- hmmStateVars(obs=obs, alpha=fvar[[1]], beta=bvar, p=seqp, b=lambda[[3]])

		# estimate model parameters - only need to estimate initital distribution and transition matrix
		lambda[est] <- hmmEstParam(gamma=svar[[1]], psi=svar[[2]], estimate=est, obs=obs)[est]

		# calculate log likelihood of the data
		logPX <- hmmProbObs(c=fvar[[2]])

		# check for terminating conditions
		if (iter > 0) {
			pdiff <- logPX - oldlogPX
			cat("Likelihood [", iter, "]: ", logPX, ", Difference in log probability: ", pdiff, "\n", sep='')
			if (pdiff < pdifflimit) break;
		} else cat("Likelihood[0]: ", logPX, "\n")

		# prepare for next iteration
		oldlogPX <- logPX
		iter <- iter + 1
	}
	
	return(lambda)
}

hmmViterbi <- function(pi, p, b, obs, steps) {
	# Viterbi algorithm
	# log scaling for underflow protection
	
	# pi: intitial probability vector
	# p: transition probability matrix
	# b: emission probability matrix
	# obs: observations
	
	nstates <- length(pi)
	T <- length(obs[[1]])
	
	# take log of model parameters
	logpi <- log(pi)
	logp <- seqPMatrix(p=p, sites=steps, logscale=1)
	logb <- list(log(b[[1]]), log(b[[2]]))
	
	# initialize data structures for storing viterbi variables
	v <- matrix(data=NA, nrow=nstates, ncol=T) # path probability matrix
	backptr <- matrix(data=NA, nrow=nstates, ncol=T) # traceback pointer matrix
	q <- rep(NA, T) # vector of states
	
	# initialization
	for (i in 1:nstates) {
		v[i,1] <- logpi[i] + logb[[i]][ obs[[1]][1], obs[[2]][1] ]
		backptr[i,1] <- 0
	}
	
	# recursion for v_t(j) = max_{q1...qt-1}(q1,...,qt-1, qt=Sj, O1,...,Ot)
	for (t in 2:T) {
		for (j in 1:nstates) {
			
			# find max(v_t(i), pij) <- should be a seperate function in C implementation
			statemax <- -Inf
			for (i in 1:nstates) {
				a <- v[i,t-1] + logp[[t]][i, j]
				if  (a > statemax) {
					statemax <- a
					backptr[j, t] <- i
				}
			}
			
			# calculate v_t(j)
			v[j, t] <- logb[[j]][ obs[[1]][t], obs[[2]][t] ] + statemax
		}
	}
	
	pstates <- max(v[,T])[1] # P(O,Q), index in case of tie
	
	# traceback to recover states that maximize P(Q|O)
	q[T] <- which(v[,T] == max(v[,T]))[1] # index in case of tie
	
	for(t in T:2) q[t-1] <- backptr[q[t],t]
	
	# print info
	cat("logP(observations, states) = ", pstates, "\n", sep='')
	
	return(q)
}

dupCoordinates <- function (q, sites) {
	# output duplicated regions, 1-indexed
	
	# q: vector of states that maximizes P(Q|O)
	# sites: matrix of sequence IDs and positions corresponding to q
	
	#regions <- data.frame(id=rep(NA,length(q)), start=rep(NA,length(q)), end=rep(NA,length(q)))
	#states <- data.frame(id=sites$V1, pos=sites$V2, state=rep(NA,length(q))) # for debugging
	#
	#states[1,3] <- ifelse(q[1] == 2, 'DUP', 'ND') # debug only
	#
	#r <- 1
	#if (q[1] == 2) {
	#	regions[1,1] <- as.character(sites[1,1])
	#	regions[1,2] <- sites[1,2]
	#}
	#
	#for (i in 2:length(q)) {
	#	states[i,3] <- ifelse(q[i] == 2, 'DUP', 'ND') # debug only
	#	
	#	# nonduplicated -> duplicated
	#	if (q[i-1] == 1 && q[i] == 2) {
	#		regions[r,1] <- as.character(sites[i,1])
	#		regions[r,2] <- sites[i,2]
	#		next
	#	}
	#
	#	# duplicated -> nonduplicated
	#	if (q[i-1] == 2 && q[i] == 1) {
	#		regions[r,3] <- sites[i-1, 2]
	#		r <- r+1
	#	}
	#}
	#
	#regions <- na.omit(regions)
	
	# faster in R:
	
	lag <- c(tail(q,-1), NA)
	breaks <- which((q == lag) == FALSE)
	start <- NA
	end <- NA
	
	# calculate regions
	regions <- NULL
	if (length(breaks) > 0) {
		if (q[breaks[1]] == 1) { # start is nonduplicated
			start <- breaks[c(TRUE,FALSE)]+1
			end <- breaks[c(FALSE,TRUE)]
			endlen <- length(end)
			if (endlen == 1 && is.na(end)) end <- length(q)
			if (start[length(start)] > end[endlen]) end <- c(end, length(q))
		} else { # start is duplicated
			start <- c(1, breaks[c(FALSE,TRUE)]+1)
			end <- breaks[c(TRUE,FALSE)]
			startlen <- length(start)
			if (is.na(start[startlen])) start <- start[-startlen]
			if (start[length(start)] > end[length(end)]) end <- c(end, length(q))
		}
		
		nregions <- length(start)
		regions <- matrix(NA, nrow=nregions, ncol=3)
		for (i in 1:nregions) {
			regions[i,1] <- as.character(sites[start[i],1])
			regions[i,2] <- sites[start[i],2]
			regions[i,3] <- sites[end[i],2]
		}
	} else {
		if (q[1] == 2) { # all sites duplicated
			regions <- matrix(NA, nrow=1, ncol=3)
			regions[1,1] <- as.character(sites[1,1])
			regions[1,2] <- sites[1,2]
			regions[1,3] <- sites[nrow(sites),2]
		}
	}
	
	states <- data.frame(id=sites$V1, pos=sites$V2, state=q-1) # 0 = nonduplicated, 1 = duplicated
	
	return(list(regions, states))
}

formatParams <- function(lrpar, covpar=NULL) {
	params <- c(lrpar, covpar)
	parmatrix <- matrix(params, nrow=1, ncol=length(params), dimnames=list(NULL,names(params)))
	return(parmatrix)
}

mainDupHmm <- function (lr, coverage=NULL, emissiontype, maxiter=100, probdiff=1e-4, lrquantile=1.0, maxcoverage=NULL, altcovmin=0, 
		lrpenalty=NULL, sampsize=NULL, paramsOnly=FALSE, emitparams=NULL) {
	# main function for duplication HMM

	# lr: ngsParalog calcLR likelihood ratios output
	# coverage: vector of average individual coverage for each site in lr
	# emissiontype: use (0) only LR or (2) LR and coverage as emissions
	# maxiter: maximum number of iterations for Baum-Welch EM
	# probdiff: minimum difference in log likelihoods between iterations of Baum-Welch EM
	# lrquantile: ignore LR values above 'lrquantile' when fitting alternative LR distribution (avoid influencce of many very extreme values)
	# maxcoverage: maximum coverage when fitting coverage distribution
	# altcovmin: lower bound for duplicated sites truncated normal coverage distribution
	# lrpenalty: type of penalty to apply to the likelihood ratio
	# sampsize: sample size (number of diploid individuals)

	rv <- list(NULL, NULL, NULL, NULL, NULL) # used to store results
	
	# penalize likelihood ratios
	
	if (! is.null(lrpenalty)) lr$V5 <- penalizeLR(lr=lr$V5, penalty=lrpenalty, sampn=sampsize)
	
	# fit emission density parameters
	
	obsfit <- list(NULL, NULL)
	if (is.null(emitparams)) {
		obsfit <- fitEmissions(lr=lr$V5, coverage=coverage, emittype=emissiontype, lrmax_quantile=lrquantile, maxcoverage=maxcoverage, min_alt_cov=altcovmin)
	} else {
		obsfit[[1]] <- emitparams[1:2]
		if (emissiontype) obsfit[[2]] <- emitparams[3:8]
	}
	
	rv[[3]] <- formatParams(lrpar=obsfit[[1]], covpar=obsfit[[2]]) # emission density parameter estimates
	if (paramsOnly) return(rv)
	
	# initialize model parameter guesses for Baum-Welch estimation
	
	lambda <- list(NA, NA, NA)
	lambda[[1]] <- initializePi() # initial distribution vector
	lambda[[2]] <- initializeP() # transitition probability matrix
	emit <- initializeEmissions(lr=lr$V5, lrpar=obsfit[[1]], coverage=coverage, covpar=obsfit[[2]], emittype=emissiontype)
	if (is.character(emit[[1]])) stop(emit[[1]],call.=FALSE) # implement better exception handling here
	lambda[[3]] <- emit[[1]] # emission probability matrix
	
	# estimate hmm parameters with Baum-Welch
	cat("\nestimating hmm parameters\n")
	estparam <- c(1,2) # only estimate initial distribution and transition probabilities
	lambda[estparam] <- hmmBaumWelch(obs=emit[[2]], steps=lr$V2, pi=lambda[[1]], p=lambda[[2]], b=lambda[[3]], maxiter=maxiter, pdifflimit=probdiff, est=estparam)[estparam]
	
	# perform decoding with Viterbi
	cat("\nfinding state sequence\n")
	q <- hmmViterbi(pi=lambda[[1]], p=lambda[[2]], b=lambda[[3]], obs=emit[[2]], steps=lr$V2)
	
	# find coordinates of duplicated regions
	cat("\ncalculating regions\n")
	regions <- dupCoordinates(q=q, sites=lr[,1:2])
	
	# organize items to return

	rv[[1]] <- regions[[1]] # duplicated regions
	rv[[2]] <- regions[[2]] # SNP states
#	rv[[4]] <- lambda[[1]] # initital state distribution - debug
#	rv[[5]] <- lambda[[2]] # transition matrix - debug
	
	return(rv)
}

parseInput <- function(input) {
	# parse parameters
	emission_type <- as.numeric(input$emit)
	if (! emission_type %in% c(0,1)) stop(paste("--emit", emission_type, "invalid. Must be 0 (LR only) or 1 (LR and coverage)"))
	
	em_max_iter <- as.numeric(input$maxiter)
	if (em_max_iter < 1) stop("--maxiter must be >= 1")
	
	pdiff <- as.numeric(input$probdiff)
	if (pdiff <= 0) stop("--probdiff must be > 0")
	
	lrcutoff <- as.numeric(input$lrquantile)
	if (lrcutoff <= 0 || lrcutoff > 1) stop("--lrquantile out of range (0,1]")
	
	maxcov <- NULL
	if (! is.null(input$maxcoverage)) {
		maxcov <- as.numeric(input$maxcoverage)
		if (maxcov <= 0) stop("--maxcoverage must be > 0")
	}
		
	alt_cov_lb <- as.numeric(input$dupcovmin)
	if (alt_cov_lb < 0) stop("--dupcovmin must be >= 0")
	
	if (!is.null(maxcov) && !is.null(alt_cov_lb) && alt_cov_lb > maxcov) stop("--dupcovmin must be <= --maxcoverage")

	samplesz <- NULL
	if (! is.null(input$n)) {
		samplesz <- ceiling(as.numeric(input$n))
		if (samplesz < 1) stop("--n must specify a sample size >= 1")
	}
	
	lrpenalty <- NULL
	if (input$penalty %in% c('none', 'NONE')) {
		lrpenalty <- NULL
	} else if (input$penalty %in% c('aic', 'AIC')) {
		lrpenalty <- 'aic'
	} else if (input$penalty %in% c('bic', 'BIC')) {
		lrpenalty <- 'bic'
		if (is.null(input$n)) stop("Use of BIC requires that the diploid sample size be provided with --n")
	} else stop(paste("--penalty", input$penalty, "invalid: Must be aic, bic, or none"))
	
	paronly <- as.numeric(input$paramOnly)
	if (! paronly %in% c(0,1)) stop(paste("--paramOnly", paronly, "invalid. Must be 1 for only estimating parameters or 0"))
	
	printstate <- as.numeric(input$printStates)
	if (! printstate %in% c(0,1)) stop(paste("--printStates", printstate, "invalid. Must be 1 for printing SNP states or 0"))
	
	# open file connections
	lrf <- read.table(input$lrfile)
	
	cov <- NULL
	if (!is.null(input$covfile)) {
		cov <- read.table(input$covfile)$V3
		if (nrow(lrf) != length(cov)) stop("Number of sites in likelihood ratio and coverage files differ")
	}
	
	emitpar <- NULL
	if (!is.null(input$paramfile)) {
		partmp <- read.table(input$paramfile, head=TRUE)
		emitpar <- as.numeric(partmp[1,])
		names(emitpar) <- colnames(partmp)
	}
	
	# perform some argument checks
	if (!is.null(emitpar)) {
		if (emission_type == 0 && length(emitpar) < 2) stop("Missing parameters in paramfile")
		if (emission_type == 1 && length(emitpar) < 8) stop("Not enough parameters supplied to use --emit 1; check paramfile")
	}
	
	rf_name <- paste(input$outfile, ".rf", sep='')
	sf_name <- paste(input$outfile, ".states", sep='')
	parf_name <- paste(input$outfile, ".par", sep='')
	
	parf_out <- NULL
	if (is.null(input$paramfile)) parf_out <- file(parf_name, "w")

	rf_out <- NULL
	if (!paronly) rf_out <- file(rf_name, "w")
	
	sf_out <- NULL
	if (!paronly && printstate) sf_out <- file(sf_name, "w")
	
	return(list(lrf, cov, rf_out, sf_out, emission_type, em_max_iter, pdiff, lrcutoff, maxcov, alt_cov_lb, lrpenalty, samplesz, paronly, emitpar, parf_out))
}

###### end functions ######

v <- paste('dupHMM.R 0.5.1',"\n") # version 4/1/2018

# parse arguments

'
Description:
Infer regions of duplication from ngsParalog likelihood ratios and sequencing depth
		
Usage:
   dupHMM.R --lrfile=<likelihood ratio file> --outfile=<output prefix> [options]
   dupHMM.R -h | --help
   dupHMM.R --version
		
Options:
   --emit=<0|1>           Use (0) only LRs or (1) LRs and coverage as emissions [default: 1]            
   --covfile=<file>       File with the average individual coverage for all sites in the LR file
   --penalty=<character>  Penalize likelihood ratios using AIC, BIC, or none for no penalty [default: none]
   --n=<int>              Diploid sample size (required for BIC)
   --maxiter=<int>        Maximum number of Baum-Welch iterations [default: 100]
   --probdiff=<float>     Minimum difference in log likelihood between Baum-Welch iterations [default: 1e-4]
   --lrquantile=<float>   Ignore LRs above this quantile when fitting alternative LR distribution [default: 1.0]
   --maxcoverage=<float>  Maximum coverage for fitting coverage distribution
   --dupcovmin=<float>    Lower bound for duplicated coverage distribution [default: 0]
   --paramOnly=<0|1>      Only estimate emission density params if 1 [default: 0]
   --paramfile=<file>     File with emission density parameter estimates
   --printStates=<0|1>    Output SNP states if 1 [default: 0]
' -> doc

opts <- docopt(doc, version=v)

# parse input
userin <- parseInput(opts)

# run the hmm
result <- mainDupHmm(lr=userin[[1]], coverage=userin[[2]], emissiontype=userin[[5]], maxiter=userin[[6]], probdiff=userin[[7]], lrquantile=userin[[8]], 
		maxcoverage=userin[[9]], altcovmin=userin[[10]], lrpenalty=userin[[11]], sampsize=userin[[12]], paramsOnly=userin[[13]], emitparams=userin[[14]])

# output results

cat("\nWriting results\n")

# regions
if (! is.null(userin[[3]])) {
	write.table(result[[1]], file=userin[[3]], row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
	close(userin[[3]])
}

# SNP states
if (! is.null(userin[[4]])) {
	write.table(result[[2]], file=userin[[4]], row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
	close(userin[[4]])
}

# emission density parameters
if (! is.null(userin[[15]])) {
	write.table(result[[3]], file=userin[[15]], row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
	close(userin[[15]])
}

cat("Finished\n")
