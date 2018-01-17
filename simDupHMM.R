# simDupHMM.R

library(expm)

simDupHmm <- function (nsnp, theta, nullr, duplr, nulcov, dupcov, pi, p, seqid='chr1') {
	# simulate paralogous data according to HMM parameters
	
	# nsnp: number of SNPs to simulate
	# theta: population scaled mutation rate (theta = 4*N*mu)
	# nullr: ngsParalog output for non-paralogous sites (should have at least nsnps sites)
	# duplr: ngsParalog output for duplicated sites (should have at least nsnps sites)
	# nulcov: matrix of total site coverage corresponding to nullr sites
	# dupcov: matrix of total site coverage corresponding to dupulr sites
	# pi: initial distribution vector [P(non-duplicated), P(duplicated)]
	# p: transition probability matrix: 
	# [nondup, nondup][nondup, dup]
	# [dup,    nondup][dup,    dup]
	
	# seqid: sequencing ID

	simsites <- data.frame(seqid = rep(NA, nsnp), site=rep(NA, nsnp), nullike=rep(NA, nsnp), altlike=rep(NA, nsnp), likeratio = rep(NA, nsnp))
	states <- data.frame(seqid = rep(NA, nsnp), pos=rep(NA, nsnp), state=rep(NA, nsnp)) # states: ND = not duplicated, DUP = duplicated
	coverage <- data.frame(seqid = rep(NA, nsnp), site=rep(NA, nsnp), coverage=rep(NA, nsnp))
	
	s <- NA # state, 1=nonduplicated, 2=duplicated
	rand <- runif(n=nsnp)
	site <- 1
	
	# find first SNP position
	while (runif(1) > theta) site <- site + 1
	
	# find initial state at first SNP
	if (rand[1] < pi[1]) {
		simsites[1,] <- c(seqid, site, nullr[1,3], nullr[1,4], nullr[1,5])
		states[1,] <- c(seqid, site, 0)
		coverage[1,] <- c(seqid, site, nulcov[1,3])
		s <- 1
	} else {
		simsites[1,] <- c(seqid, site, duplr[1,3], duplr[1,4], duplr[1,5])
		states[1,] <- c(seqid, site, 1)
		coverage[1,] <- c(seqid, site, dupcov[1,3])
		s <- 2
	}
	
	# simulate 2 .. number snps
	for (t in 2:nsnp) {
		# find SNP position
		repeat {
			site <- site + 1;
			if (runif(1) < theta) break;
		}
		
		# calc probability of being in state S_t
		nstepP <- p %^% (site - as.numeric(simsites[t-1,2]))
		
		# determine state
		if (s == 1) {
			if (rand[t] < nstepP[s,2]) s <- 2
		} else {
			if (rand[t] < nstepP[s,1]) s <- 1
		}
		
		# draw LR for the current state
		if (s == 1) {
			nn <- sample(1:nrow(nullr), 1)
			simsites[t,] <- c(seqid, site, nullr[nn,3], nullr[nn,4], nullr[nn,5])
			states[t,] <- c(seqid, site, 0)
			coverage[t,] <- c(seqid, site, nulcov[nn,3])
		} else {
			dn <- sample(1:nrow(duplr), 1)
			simsites[t,] <- c(seqid, site, duplr[dn,3], duplr[dn,4], duplr[dn,5])
			states[t,] <- c(seqid, site, 1)
			coverage[t,] <- c(seqid, site, dupcov[dn,3])
		}
	}
	
	# return simulated LRs and states
	return(list(simsites, states, coverage))
}

subBalancingSnps <- function (subrange, simsites, states, coverage, ballr, balcov) {
	# changes nonduplicated sites in 'range' to nonduplicated balancing selection sites

	# subrange: vector of start and stop positions of balanacing selection sites
	# simsites: dataframe of duplication likelihood ratios, this should be simDupHmm[[1]]
	# states: dataframe of HMM states corresponding to simsites, this should be simDupHmm[[2]]
	# coverage: dataframe of coverage corresponding to simsites, this should be simDupHmm[[3]]
	# ballr: ngsParalog output of balancing selection sites to be substitutded in
	# balcov: dataframe of total site coverage corresponding to ballr

	if (length(subrange) %% 2 != 0) stop("Odd number of stop/start coordinates for range")
	
	idx <- rep(NA, length(subrange))
	for (i in 1:length(subrange)) idx[i] <- which(states$pos == subrange[i])

	# substitute balancing selection sites for nonduplicated sites
	for (i in seq(from=1, to=(length(idx)-1), by=2)) {
		# make sure substitution region does not overlap a duplicated region
		if (length(which(states$state[idx[i]:idx[i+1]]==1)) > 0) stop(paste(subrange[i], "-", subrange[i+1], " spans a duplicated region", sep=''))
		
		# substitute
		rand <- sample(x=1:nrow(ballr), size=(idx[i+1]-idx[i])+1)
		k <- 1;
		for (j in idx[i]:idx[i+1]) {
			simsites[j,3:5] <- ballr[rand[k], 3:5];
			states[j,3] <- 2 # let 2 = "BAL"
			coverage[j,3] <- balcov[rand[k],3] 
			k <- k+1
		}
	}

	return(list(simsites, states, coverage))	
}

hmmError <- function(simstates, hmmstates) {
	# simstates: simulated states
	# hmmstates: states inferred with dupHMM.R

	simstates <- as.character(simstates)
	hmmstates <- as.character(hmmstates)
	simstates <- replace(simstates, simstates == 2, 0)
	n <- length(simstates)
	falsepos <- 0
	falseneg <- 0
	for (i in 1:n) {
		if (simstates[i] != hmmstates[i]) {
			if (simstates[i] == 0) falsepos <- falsepos+1 else falseneg <- falseneg+1		
		}
	}
	misid <- falsepos + falseneg
	error <- c(misid/n, falsepos/n, falseneg/n)
	names(error) <- c("mis_rate", "false_positive", "false_negative")
	return(error) 
}
