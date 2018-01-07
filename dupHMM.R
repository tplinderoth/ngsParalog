# dupHMM.R

# TODO:
###############################################################################

library(expm, quietly=TRUE)
library(truncnorm, quietly=TRUE)
library(argparser, quietly=TRUE)

###### functions ######

initializePi <- function() {
	# initialize initial distribution vector uniformly
	
	pi <- c(0.5, 0.5)
	return(pi)
}

initializeP <- function() {
	# initialize transition probability matrix uniformly
	
	nstates = 2
	trans <- matrix(NA, nrow=nstates, ncol=nstates)

	p <- 1/nstates
	for (i in 1:nstates) {
		for (j in 1:nstates) {
			trans[i,j] <- p
		}
	}
	
	return(trans)
}

fitlr <- function (lr, coverage, nullmean, nullsd, tailcutoff = 1.0) {
	# lr: likelihood ratios
	# coverage: average individual coverage
	# nullmean: estimate of average individual coverage for nonduplicated sites
	# nullsd: estimate of the standard deviation for the average individual coverage for nonduplicated sites
	# tailcutoff: fit alternative distribution to LRs under the 'tailcutoff' quantile (avoid influence of extreme values)
	
	# subset of LR values for sites that have < 0.95 probability of coming from null according to coverage and are within the 'tailcutoff' percentile
	covprob <- 0.95 # coverage used to avoid influencing fit by null sites with heterozygote advantage
	covlower <- qtruncnorm(covprob, mean=nullmean, sd=nullsd, a=0, b=Inf)
	lrupper <- quantile(lr, tailcutoff)
	sublr <- lr[which(coverage > covlower & lr < lrupper)]
	
	# estimate proportion of LRs under the null
	lrquantile <- 2.705 # 95% quantile for null LR distribution
	pnull <- 1 - length(which(sublr > lrquantile))/length(sublr)
	
	# approximate median of alternative LR values
	ncpguess <- mean(sublr)
	
	fit <- optim(par=c(pnull, 1, ncpguess), fn=lrllh, method="L-BFGS-B", lower=c(0, 0, 0), upper=c(1, Inf, Inf), x=sublr)
	
	return(fit$par)
}

lrllh <- function (par, x) {
	# LR ~ par[1]*[0.5*chisq(df=0, ncp=0) + 0.5*chisq(df=1, ncp=0)] + (1-par[1])*chisq(df=par[2], ncp=par[3])
	# par[1] = probability LR comes from null distribution 
	# par[2] = df for alternate chisquare distribution
	# par[3] = noncentrality parameter for alternate chisquare distribution
	
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
	y2 <- dzero*dchisq(x, df=1, ncp=0) + (1-par[1])*dchisq(x, df=par[2], ncp=par[3])
	
	llh <- c(y1, y2)
	llh <- replace(llh, llh==0, .Machine$double.xmin)
	-sum(log(llh))
}

covllh <- function (par, x) {
	# null coverage ~ truncated_normal(mean=par[1], sd=par[2], lowerbound=0, upperbound=Inf)
	# par[1] = truncated normal mean
	# par[2] = truncated normal standard deviation
	
	llh <- dtruncnorm(x, mean=par[1], sd=par[2], a=0, b=Inf)
	llh <- replace(llh, llh==0, .Machine$double.xmin)
	-sum(log(llh))
}

fitCoverage <- function(coverage, lr) {
	# approximate mean and variance of average individual coverage
	nullcutoff <- 1
	nondupidx <- which(lr<nullcutoff)
	nullavg <- mean(coverage[nondupidx])
	nullsd <- sqrt(var(coverage[nondupidx]))
	
	# approximate proportion of nonduplicated sites
	pnull <- length(nondupidx)/length(lr)
	
	# parameter bounds
	mu_min <- 0
	mu_max <- 3*nullavg
	sd_min <- 0
	sd_max <- Inf
	
	# find MLE parameter values
	fit <- optim(par=c(nullavg, nullsd), fn=covllh, method="L-BFGS-B", lower=c(mu_min, sd_min), upper=c(mu_max, sd_max), x=coverage)
	
	return(fit$par)
}

initializeEmissions <- function(lr, coverage, lrmax_quantile) {
	# initialize emission probabilites
	
	# lr: vector of duplicaiton likelihood ratios
	# coverage: vector of average individual coverage
	
	lrseq <- seq(from=2, to=max(ceiling(lr)), by=1)
	covseq <- seq(from=2, to=max(ceiling(coverage)), by=1)
	
	lrprob <- matrix(NA, nrow=2, ncol=length(lrseq)+2)
	covprob <- matrix(NA, nrow=2, ncol=length(covseq)+2)

	# estimate parameters for null coverage distribution
	cat("Fitting coverage distribution\n")
	covpar <- fitCoverage(coverage, lr)
	
	# estimate parameters for LR distribution
	cat("Fitting LR distribution\n")
	lrpar <- fitlr(lr=lr, coverage=coverage, nullmean=covpar[1], nullsd=covpar[2], tailcutoff=lrmax_quantile)

	# calculate emission probabilities
	cat("Calculating emission probabilities\n")
	
	# calculate probabilites of the observed LR
	lrmax <- 1
	seenalt <- 0
	lrprob[1,1] <- 0.5 + 0.5*pchisq(1,1)
	lrprob[2,1] <- pchisq(1, df=lrpar[2], ncp=lrpar[3])
	nullsum <- lrprob[1,1]
	altsum <- lrprob[2,1]
	
	for (i in lrseq) {
		lrprob[1,i] <- (0.5 + 0.5*pchisq(q=i, df=1)) - nullsum # null
		lrprob[2,i] <- pchisq(i, df=lrpar[2], ncp=lrpar[3]) - altsum # alternate
		nullsum <- nullsum + lrprob[1,i]
		altsum <- altsum + lrprob[2,i]
		if (seenalt < 1 && lrprob[2,i] > 0) seenalt <- 1 # check if the alternative distribution has been entered 
		if (lrprob[1,i] == 0 && lrprob[2,i] == 0 && seenalt > 0) break
		lrmax <- lrmax + 1
	}
	
	# assign very small probability to LR emissions with zero probability under the null and alternative distributions
	lrprob[,which(lrprob[1,]==0 & lrprob[2,]==0)] <- .Machine$double.xmin
	
	# calculate probabilities of the observed coverage
	covmax <- 1
	seenalt <- 0
	covmean_alt <- 2*covpar[1] # assume duplicated sites have twice the average coverage as nonduplicated sites
	covprob[1,1] <- ptruncnorm(1, mean=covpar[1], sd=covpar[2], a=0, b=Inf)
	covprob[2,1] <- ptruncnorm(1, mean=covmean_alt, sd=covpar[2], a=0, b=Inf)
	nullsum <- covprob[1,1]
	altsum <- covprob[2,1]
	
	for (i in covseq) {
		covprob[1,i] <- ptruncnorm(i, mean=covpar[1], sd=covpar[2], a=0, b=Inf) - nullsum # null
		covprob[2,i] <- ptruncnorm(i, mean=covmean_alt, sd=covpar[2], a=0, b=Inf) - altsum # alternate
		nullsum <- nullsum + covprob[1,i]
		altsum <- altsum + covprob[2,i]
		if (seenalt > 0 && covprob[2,i] > 0) seenalt <- 1 # check if alternative distribution has been entered
		if (covprob[1,i] == 0 && covprob[2,i] == 0 && seenalt > 0) break
		covmax <- covmax + 1
	}
	
	# assign very small probability to coverage emissions with zero probability under the null and alternative distributions
	covprob[,which(covprob[1,]==0 & covprob[2,]==0)] <- .Machine$double.xmin
	
	# calculate joint probabilities of LR and coverage
	b <- list(matrix(NA, nrow=lrmax, ncol=covmax), matrix(NA, nrow=lrmax, ncol=covmax))
	for (i in 1:lrmax) {
		for (j in 1:covmax) {
			b[[1]][i,j] <- lrprob[1,i] * covprob[1,j] # null
			b[[2]][i,j] <- lrprob[2,i] * covprob[2,j] # alternative
		}
	}
	
	# make the emission probs sum to one
	b[[1]] <- b[[1]]/sum(b[[1]]) # null emissions must sum to 1
	b[[2]] <- b[[2]]/sum(b[[2]]) # alternative emissions must sum to 1
	
	# bound the maximum LR and coverage values
	lr[lr==0] <- 0.1
	lr <- ceiling(lr)
	lr <- replace(lr, lr > lrmax, lrmax)
	
	coverage[coverage==0] <- 0.1
	coverage <- ceiling(coverage)
	coverage <- replace(coverage, coverage > covmax, covmax)
	
	return(list(b, list(lr, coverage)))
}

seqPMatrix <- function(p, sites, logscale=0) {
	# calculates a list of n-step transition matrices based on the distance between snps

	# p: transition probability matrix
	# sites: vector steps of the HMM (SNP positions)
	
	P <- list(p) # P[[1]] is a dummy
	if (logscale) P[[1]] = log(P[[1]])
	
	for (i in 2:length(sites)) {
		P[[i]] = p %^% (sites[i] - sites[i-1])
		if (logscale) P[[i]] = log(P[[i]])
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
	
	nstates = length(pi)
	T = length(obs[[1]])
	alpha <- matrix(data=NA, nrow=T, ncol=nstates) # alpha(t,i) = P(O1, O2, O3,..., Ot, qt=Si), forward variables
	c <- rep(NA, T) # scaling variables
	
	# initialization step - might need to adjust this and start at alpha(0,i) as in Durbin
	c[1] = 0
	for (i in 1:nstates) {
		alpha[1, i] = pi[i] * b[[i]][ obs[[1]][1], obs[[2]][1] ] # alpha(1,i) = pi_i * b_i(O1)
		c[1] = c[1] + alpha[1,i]
	}
	
	# scale the alpha[1,i]
	for (i in 1:nstates) alpha[1, i] = alpha[1, i]/c[1]

	# induction step
	for (t in 2:T) {
		c[t] = 0
		
		# alpha(t+1, i) = sum_i=1_to_N( alpha(t,j) * p_ji * b_i(O_t+1) )
		for (i in 1:nstates) {
			alpha[t, i] = 0
			for (j in 1:nstates) {
				alpha[t, i] = alpha[t, i] + alpha[t-1, j] * p[[t]][j, i]
			}
			alpha[t, i] = alpha[t, i] * b[[i]][ obs[[1]][t], obs[[2]][t] ]
			c[t] = c[t] + alpha[t, i]
		}
		
		# scale the forward variables
		for (i in 1:nstates) alpha[t, i] = alpha[t, i]/c[t]
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
	
	nstates = nrow(p[[1]])
	T = length(obs[[1]])
	beta <- matrix(data=NA, nrow=T, ncol=nstates) # beta(t,i) = P(O_t+1, O_t+2,..., O_T | q_t = S_i), backward variables
	
	# initializaiton
	for (i in 1:nstates) {
		beta[T, i] = 1/c[T]
	}
	
	# recursion
	for (t in (T-1):1) {
		for (i in 1:nstates) {
			beta[t, i] = 0
			for (j in 1:nstates) {
				beta[t, i] = beta[t, i] + p[[t+1]][i, j] * b[[j]][ obs[[1]][t+1], obs[[2]][t+1] ] * beta[t+1, j]
			}
			beta[t, i] = beta[t, i]/c[t]
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
	gamma <- lapply(1:T-2, matrix, data=NA, nrow=nstates, ncol=nstates) # P(x_t = i, x_t+1 = j | x, lambda)
	psi <- matrix(0, nrow=T, ncol=nstates) # P(x_t = i | x, lambda)
	
	for (t in 1:(T-1)) {
		
		# calculate denominator, P(X)
		px = 0
		for (i in 1:nstates) {
			for (j in 1:nstates) px = px + ( alpha[t,i] * p[[t+1]][i, j] * b[[j]][ obs[[1]][t+1], obs[[2]][t+1] ] * beta[t+1, j] )
		}
		
		
		for (i in 1:nstates) {
			psi[t, i] = 0
			for (j in 1:nstates) {
				# probability that p_ij is used at position t
				gamma[[t]][i,j] = (alpha[t, i] * p[[t+1]][i, j] * b[[j]][ obs[[1]][t+1], obs[[2]][t+1] ] * beta[t+1, j]) / px
				# probability that x_t = i
				psi[t, i] = psi[t, i] + gamma[[t]][i,j]
			}
		}
	}
	
	px = 0
	for (i in 1:nstates)
		px = px + alpha[T, i]
	for (i in 1:nstates)
		psi[T, i] = alpha[T, i]/px
	
	return(list(gamma, psi))
}

hmmEstParam <- function(gamma, psi, estimate=c(1,2,3), obs=NA) {
	# Baum-Welch parameter re-estimation step 
	# obs required if estimating emission matrix
	
	# gamma: probility that p_ij is used at step t
	# psi: P(x_t) = i at each step t
	# estimate: what parameters to restimate, 1=>pi, 2=>transition probabilities, 3=>emission probabilities
	# obs: vector of potential observations

	nstates <- ncol(psi)
	T <- nrow(psi)
	params <- list(NA, NA, NA) # [[1]] = pi, [[2]] = transition matrix, [[3]] = emission matrix
	
	# estimate initial distribution
	if (1 %in% estimate) {
		params[[1]] <- rep(NA, nstates)
		for (i in 1:nstates) params[[1]][i] = psi[1, i]
	}
	
	# estimate transition probability matrix
	if (2 %in% estimate) {
		params[[2]] <- matrix(NA, nrow=nstates, ncol=nstates)
		
		for (i in 1:nstates) {
			for (j in 1:nstates) {
				expij = 0 # expected number of transitions from state i to j
				exi = 0 # expected number of transitions from state i - can store this for estimating B
				for (t in 1:(T-1)) {
					expij = expij + gamma[[t]][i,j]
					exi = exi + psi[t,i]
				}
				params[[2]][i,j] = expij/exi
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
			exi = 0 # expected number of times in state i
			for (t in 1:T) {
				params[[3]][[i]][ obs[[1]][t], obs[[2]][t] ] <- params[[3]][[i]][ obs[[1]][t], obs[[2]][t] ] + psi[t, i]
				exi = exi + psi[t, i]
			}
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
	
	pobs = 0
	for (j in 1:length(c)) pobs = pobs + log(c[j])
	
	return(pobs)
}

hmmBaumWelch <- function(obs, steps, pi, p, b, maxiter=100, pdifflimit=1e-4, est=c(1,2,3)) {
	# Baum Welch algorithm
	# only returns updated estimates of pi and p
	
	# obs: observations
	# pi: initital probability distribution
	# p: transition probability matrix
	# b: emission probability matrix
	# maxiter: maximum number of EM iterations
	# logdiff: minimum difference in the log likelihoods
	# parameters to estimate: 1=>inititial distribution, 2=>transition matrix, 3=>emission probs
	
	lambda = list(pi, p, b)
	oldlogPX = 0
	logPX = 0
	iter = 0
	
	while (iter < maxiter) {
		seqp = seqPMatrix(p=lambda[[2]], sites=steps, logscale=0)
		
		# calculate forward variables
		fvar = hmmForward(obs=obs, pi=lambda[[1]], p=seqp, b=lambda[[3]])
		
		# calculate backward variables
		bvar = hmmBackward(obs=obs, p=seqp, b=lambda[[3]], c=fvar[[2]])
		
		# calculate probability of i->j and of being in state i at time t
		svar = hmmStateVars(obs=obs, alpha=fvar[[1]], beta=bvar, p=seqp, b=lambda[[3]])
		
		# estimate model parameters - only need to estimate initital distribution and transition matrix
		lambda[est] = hmmEstParam(gamma=svar[[1]], psi=svar[[2]], estimate=est, obs=obs)[est]
		
		# calculate log likelihood of the data
		logPX = hmmProbObs(c=fvar[[2]])
		
		#print(lambda[est]) # debug
		
		# check for terminating conditions
		if (iter > 0) {
			pdiff = logPX - oldlogPX
			cat("Likelihood [", iter, "]: ", logPX, ", Difference in log probability: ", pdiff, "\n", sep='')
			if (pdiff < pdifflimit) break;
		} else cat("Likelihood[0]: ", logPX, "\n")
		
		# prepare for next iteration
		oldlogPX = logPX
		iter = iter + 1
	}
	
	#return(lambda)
	return(lambda)
}

hmmViterbi <- function(pi, p, b, obs, steps) {
	# Viterbi algorithm
	# log scaling for underflow protection
	
	# pi: intitial probability vector
	# p: transition probability matrix
	# b: emission probability matrix
	# obs: observations
	
	nstates = length(pi)
	T = length(obs[[1]])
	
	# take log of model parameters
	logpi <- log(pi)
	logp <- seqPMatrix(p=p, sites=steps, logscale=1)
	logb <- list(log(b[[1]]), log(b[[2]])) # might have to address log(0)
	
	# initialize data structures for storing viterbi variables
	v <- matrix(data=NA, nrow=nstates, ncol=T) # path probability matrix
	backptr <- matrix(data=NA, nrow=nstates, ncol=T) # traceback pointer matrix
	q <- rep(NA, T) # vector of states
	
	# initialization
	for (i in 1:nstates) {
		v[i,1] = logpi[i] + logb[[i]][ obs[[1]][1], obs[[2]][1] ]
		backptr[i,1] = 0
	}
	
	# recursion for v_t(j) = max_{q1...qt-1}(q1,...,qt-1, qt=Sj, O1,...,Ot)
	for (t in 2:T) {
		for (j in 1:nstates) {
			
			# find max(v_t(i), pij) <- should be a seperate function in C implementation
			statemax = -.Machine$double.xmax
			for (i in 1:nstates) {
				a = v[i,t-1] + logp[[t]][i, j]
				if  (a > statemax) {
					statemax = a
					backptr[j, t] = i
				}
			}
			
			# calculate v_t(j)
			v[j, t] = logb[[j]][ obs[[1]][t], obs[[2]][t] ] + statemax
		}
	}
	
	pstates = max(v[,T])[1] # P(O,Q), index in case of tie
	
	# traceback to recover states that maximize P(Q|O)
	q[T] = which(v[,T] == max(v[,T]))[1] # index in case of tie
	
	for(t in T:2) q[t-1] = backptr[q[t],t]
	
	# print info
	cat("-logP(observations, states) = ", pstates, "\n", sep='')
	
	return(q)
}

dupCoordinates <- function (q, sites) {
	# output duplicated regions, 1-indexed
	
	# q: vector of states that maximizes P(Q|O)
	# sites: matrix of sequence IDs and positions corresponding to q
	
	regions <- data.frame(id=rep(NA,length(q)), start=rep(NA,length(q)), end=rep(NA,length(q)))
	states <- data.frame(id=sites$V1, pos=sites$V2, state=rep(NA,length(q))) # for debugging
	
	states[1,3] = ifelse(q[1] == 2, 'DUP', 'ND') # debug only
	
	r = 1
	if (q[1] == 2) {
		regions[1,1] = as.character(sites[1,1])
		regions[1,2] = sites[1,2]
	}
	
	for (i in 2:length(q)) {
		states[i,3] = ifelse(q[i] == 2, 'DUP', 'ND') # debug only
		
		# nonduplicated -> duplicated
		if (q[i-1] == 1 && q[i] == 2) {
			regions[r,1] = as.character(sites[i,1])
			regions[r,2] = sites[i,2]
			next
		}
	
		# duplicated -> nonduplicated
		if (q[i-1] == 2 && q[i] == 1) {
			regions[r,3] = sites[i-1, 2]
			r = r+1
		}
	}

	regions <- na.omit(regions)
	
	return(list(regions, states))
}

###### end functions ######

###### main ######

mainDupHmm <- function (lr, coverage, maxiter=100, probdiff=1e-4, lrquantile=0.99) {
	# main function for duplication HMM
	 
	# lr: ngsParalog calcLR likelihood ratios output
	# coverage: vector of average individual coverage for each site in lr
	# maxiter: maximum number of iterations for Baum-Welch EM
	# probdiff: minimum difference in log likelihoods between iterations of Baum-Welch EM
	# lrquantile: ignore LR values above 'lrquantile' when fitting alternative LR distribution (avoid influencce of many very extreme values) 

	# initialize model parameter guesses for Baum-Welch estimation
	
	lambda <- list(NA, NA, NA)
	lambda[[1]] <- initializePi() # initial distribution vector
	lambda[[2]] <- initializeP() # transitition probability matrix
	emit <- initializeEmissions(lr=lr$V5, coverage=coverage, lrmax_quantile=lrquantile) # returns [emission matrix, [lr, coverage]]
	lambda[[3]] <- emit[[1]] # emission probability matrix
	
	# estimate hmm parameters with Baum-Welch
	estparam <- c(1,2) # only estimate initial distribution and transition probabilities
	lambda[estparam] <- hmmBaumWelch(obs=emit[[2]], steps=lr$V2, pi=lambda[[1]], p=lambda[[2]], b=lambda[[3]], maxiter=maxiter, pdifflimit=probdiff, est=estparam)[estparam]
	
	# perform decoding with Viterbi
	q <- hmmViterbi(pi=lambda[[1]], p=lambda[[2]], b=lambda[[3]], obs=emit[[2]], steps=lr$V2)
	
	# find coordinates of duplicated regions
	regions <- dupCoordinates(q=q, sites=lr[,1:2])
	
	#return(regions) # returns a list of regions and states - only return regions in final implementation
	return(list(regions[[1]], regions[[2]], lambda[[1]], lambda[[2]])) # return states and estimate of initial state distribution and transition matrix for debugging
}

###### end main ######