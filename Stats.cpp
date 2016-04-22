/*
 * stats.cpp
 */

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <math.h>
#include <sys/stat.h>
#include <time.h>
#include "parsePileup.h"
#include "Stats.h"
#include "bfgs.h"
#include <iostream> // debug

double Stats::optimLR (Optim* null, Optim* alt, double (*fn)(const double x[], const void*), void (*dfn)(const double x[], double y[], const void*), int islog, int* status)
{
	bool weightcount = true;
	double like;
	int i = 0;
	*status = 0;

	// find MLE for f frequency under the null hypothesis

	// use the sample alternate allele frequency as a starting point
	null->par[0] = mafguess(static_cast<Pileup*>(null->data), weightcount);
	like = findmax_bfgs(null->getDim(), null->par, null, fn, dfn, null->lowbounds(), null->upbounds(), null->numbounds(), null->verblevel(), null->fail());

	// check for failure of guessed maf start point optimization
	if (null->isfail())
	{
		null->multiOptim(fn, dfn);
		// check for optimization failure at multiple start points
		if (null->isfail())
		{
			fprintf(stderr, "Null optimization failed. ");
			*status = 1;
			return 0.0;
		}
	}
	else
	{
		null->setllh() = like; //set
		for (i=0; i < null->getDim(); ++i)
			null->setmlparam(i, null->par[i]);
	}

	// find MLE for f and m under the alternative hypothesis

	// try multiple arbitrary start points
	alt->multiOptim(fn, dfn);

	// check for optimization failure at multiple start points
	if (alt->isfail())
	{
		*status = 1;
		alt->setllh() = 1.0/0.0;
		*(alt->fail()) = 0;
	}

	// use null MLE estimates as start point
	alt->par[0] = null->mlparam(0);
	alt->par[1] = 1.0;
	like = findmax_bfgs(alt->getDim(), alt->par, alt, fn, dfn, alt->lowbounds(), alt->upbounds(), alt->numbounds(), alt->verblevel(), alt->fail());

	// check for optimization failure from null MLE start point
	if (alt->isfail() && *status)
	{
		fprintf(stderr, "Alternative optimization failed. ");
		return 0.0;
	}

	if (like < alt->llh())
	{
		alt->setllh() = like;
		for (i=0; i < alt->getDim(); ++i)
			alt->setmlparam(i, null->par[i]);
	}

	// return likelihood ratio
	return( calcLR(null->llh(), alt->llh(), islog) );
}

double Stats::mafguess (Pileup* pile, bool wt)
{
	char minor;
	const static char a [] = {'A', 'C', 'G', 'T'};
	int i;
	double c = 0.0;
	double total = 0.0;
	double n;

	minor = pile->minorid() ? pile->minorid() : pile->empiricalMinorFast(wt);

	for (i=0; i<4; ++i)
	{
		n = wt ? pile->wtalleleCount(a[i]) : static_cast<double>(pile->alleleCount(a[i]));
		if (a[i] == minor) c = n;
		total += n;
	}

	return (c/total);
}


double Stats::calcLR (const double null, const double alt, bool islog)
{
//calculates likelihood ratio statistic

	double lr = 0.0;

	if (islog) // likelihoods are -log scaled
		lr = 2 * null - 2 * alt;
	else
		lr = -2 * log(null) + 2 * log(alt);
	return lr;
}

double Stats::negLogfn (const double para [], const void *generic_dat)
{
	/*
	 * main likelihood function
	 */

	const Optim* optdata = static_cast<const Optim*>(generic_dat);
	const Pileup* pile = static_cast<const Pileup*>(optdata->data);
	const int npara = 2; // number parameters in model
	const int findex = 0; // index position of allele frequency parameter
	const int mindex = 1; // index position of admixture proportion
	static double p [npara]; // para[0] = alternate allele frequency, para[1] = admixture proportion
	static double epsilon; // Phred scaled quality score
	static double geno_marg; // marginal sum over genotype2 configurations
	static double loglike; // log likelihood of likelihood function
	static double genoprior [3];
	static int k1, k2;
	static Matrix<double> indlike (3, 3); // individual llh for all genotypic configurations
	static double maxlike; // max individual llh among all genotypic configurations

	if (optdata->getDim() < npara) // null case; fix m to 1
	{
		p[findex] = para[0]; // alternate allele frequency (f)
		p[mindex] = 1.0; // admixture proportion (m)
	}
	else // alternative case
	{
		for (int i = 0; i < optdata->getDim(); ++i)
				p[i] = para[i];
	}
	loglike = 0.0;

	// identify major and minor alleles
	char major = pile->majorid();
	char minor = pile->minorid();

	// P(G1) does not depend on f: k1 can be 0 or 2
	genoprior[0] = genoPrior(p[findex], 0); // P(G1 = {0,2}, G2 = 0|f)
	genoprior[1] = genoPrior(p[findex], 1);  // P(G1 = {0,2}, G2 = 1|f)
	genoprior[2] = genoPrior(p[findex], 2); // P(G1 = {0,2}, G2 = 2|f)
	static std::vector<SiteData>::const_iterator ind_iter;
	static std::vector<seqread>::const_iterator readIter;

	for (ind_iter = pile->seqdat.begin(); ind_iter != pile->seqdat.end(); ++ind_iter) // sum over all individuals
	{
		if (ind_iter->cov() < 1) // missing data for individual
			continue;
		maxlike = -1.0/0.0;
		for (k1 = 0; k1 <= 2; k1 += 2) // sum over genotype1
		{
			for (k2 = 0; k2 <= 2; ++k2) // sum over genotype2
			{
				indlike[k1][k2] = 0.0;
				for (readIter = ind_iter->rdat.begin(); readIter != ind_iter->rdat.begin() + ind_iter->cov(); ++readIter) // product over all reads of individual
				{
					{
						epsilon = pow( 10, - (static_cast <double > (readIter->second)) / 10 );
						indlike[k1][k2] += log(readProb(p[mindex], epsilon, k1, k2, major, minor, readIter->first));
					}
				}
				if (indlike[k1][k2] > maxlike)
					maxlike = indlike[k1][k2];
			}
		}
		geno_marg = 0.0;
		for (k1 = 0; k1 <= 2; k1 += 2)
			for (k2 = 0; k2 <= 2; ++k2)
				geno_marg += exp(indlike[k1][k2] - maxlike) * genoprior[k2]; // switch back to probability space
		loglike += log(geno_marg) + maxlike;
	}
	return -loglike;
}


double Stats::readProb (const double m, const double err, const int g1, const int g2, const char ref, const char alt, const char obs)
{
	/*
	 * returns P(obs_read | allele1)*P(allele1 | G1,G2) + P(obs_read | allele2)*P(allele2 | G1,G2)
	 */

	double pread = 1.0;
	if (g1 == 0)
	{
		if (g2 == 0)
		{
			if (obs==ref)
				pread = 1.0 - err;
			else if (obs==alt)
				pread = err/3.0;
			else
				pread = err/3.0;
		}
		else if (g2 == 1)
		{
			if (obs==ref)
				pread = (1.0-err)*(1.0-0.5*m) + (err*m)/6.0;
			else if (obs==alt)
				pread = (1.0-err)*(m/2.0) + (err/3.0)*(1.0-0.5*m);
			else
				pread = (err/3.0)*(1.0-0.5*m) + (err*m)/6.0;
		}
		else if (g2 == 2)
		{
			if (obs==ref)
				pread = (1.0-err)*(1.0-m) + (err*m)/3.0;
			else if (obs==alt)
				pread = (1.0-err)*m + (err/3.0)*(1.0-m);
			else
				pread = (err/3.0)*(1.0-m) + (err*m)/3.0;
		}
		else
			genoErr(g1,g2);
	}
	else if (g1 == 2)
	{
		if (g2 == 0)
		{
			if (obs==ref)
				pread = (1.0-err)*m + (err/3.0)*(1.0-m);
			else if (obs==alt)
				pread = (1.0-err)*(1.0-m) + (err*m)/3.0;
			else
				pread = (err/3.0)*(1.0-m) + (err*m)/3.0;
		}
		else if (g2 == 1)
		{
			if (obs==ref)
				pread = (1.0-err)*(m/2.0) + (err/3.0)*(1.0-0.5*m);
			else if (obs==alt)
				pread = (1.0-err)*(1.0-0.5*m) + (err*m)/6.0;
			else
				pread = (err/3.0)*(1.0-0.5*m) + (err*m)/6.0;
		}
		else if (g2 == 2)
		{
			if (obs==ref)
				pread = err/3.0;
			else if (obs==alt)
				pread = 1.0-err;
			else
				pread = err/3.0;
		}
		else
			genoErr(g1,g2);
	}
	else
		genoErr(g1,g2);

	return pread;
}

double Stats::genoPrior (const double f, const int g2)
{
	/*
	 * calculates P(G1,G2|f)
	 */

	double prob = 0;

	if (g2 == 0)
		prob = 0.5 * (1-f) * (1-f);
	else if (g2 == 1)
		prob = f * (1-f);
	else if (g2 == 2)
		prob = 0.5 * f * f;
	else
		fprintf(stderr, "Invalid genotype2: g2 = %i", g2);

	return prob;
}

void Stats::kahanSum(double summand, double* total, double* comp)
{
	double x = summand - *comp;
	double y = *total + x;
	*comp = (y - *total) - x;
	*total = y;
}

void Stats::genoErr (const char g1, const char g2)
{
	fprintf(stderr, "Invalid genotype configuration: geno1 = %i, geno2 = %i\n", g1, g2);
}


/*
 * The following code for calculating an analytic gradient is deprecated - switched to numeric gradient.
 * It will need to be fixed and updated if ever used again.

void Stats::calcGradient (const double para [], double grad [], const void* generic_dat)
{
	// gets the gradient of the negative log likelihood function
	// para[0] = alternate allele frequency (f); para[1] = admixture proportion (M)

	const FunData* fnvals = static_cast<const FunData*>(generic_dat);
	const Pileup* pile = static_cast<const Pileup*>(fnvals->sitedat);
	static int i, j, confign, k1, k2; // counters
	const int npara = 2; // number of parameters in model
	const int findex = 0; // index position of allele frequency parameter
	const int mindex = npara - 1 - findex; // index position of admixture proportion
	static double p [npara]; // parameters
	static double correct_prec [npara]; // compensation for lost low-order bits in Kahan Sum
	static double epsilon; // Phred scaled quality score
	static double innerlogder [npara]; // [df/df, df/dm];
	static Matrix<double> geno_diff (6, npara); // holds partial derivative for each G1 G2 configuration
	static double indl; // P(individual Data| G1,G2)P(G1,G2|f)
	static double geno_marg;
	static double dread_prob [npara] = {}; // [df/df, df/dm] for sum over P(x|b)P(b|G1,G2)
	static double indsum[npara]; // holds partial derivative for individual
	static double genoprior [3];
	static double dprior [3];
	static Matrix<double> indlike (3, 3); // llh for individual at all genotypic configurations
	static double maxlike; // max llh for individual among all genotypic configurations
	static double maxprob; // maxlike in probability space

	if (fnvals->optp < fnvals->np) // null case; fix m to 1
	{
		p[findex] = para[findex]; // alternate allele frequency (f)
		p[mindex] = 1.0; // admixture proportion (m)
	}
	else // alternative case
	{
		for (int i = 0; i < npara; ++i)
			p[i] = para[i];
	}

	// P(G1) does not depend on f; k1 can be 0 or 2
	genoprior[0] = genoPrior(p[findex], 0, 0); // P(G1 = {0,2}, G2 = 0|f)
	genoprior[1] = genoPrior(p[findex], 0, 1); // P(G1 = {0,2}, G2 = 1|f)
	genoprior[2] = genoPrior(p[findex], 0, 2); // P(G1 = {0,2}, G2 = 2|f)
	dprior[0] = diffGenoPrior(p[findex], 0, 0); // d/df P(G1 = {0,2}, G2 = 0|f)
	dprior[1] = diffGenoPrior(p[findex], 0, 1); // d/df P(G1 = {0,2}, G2 = 1|f)
	dprior[2] = diffGenoPrior(p[findex], 0, 2); // d/df P(G1 = {0,2}, G2 = 2|f)

	static std::vector<SiteData>::const_iterator ind_iter;
	static std::vector< std::pair<int, int> >::const_iterator readdat;

	for (i = 0; i < fnvals->optp; ++i) // set gradient values to zero
	{
		grad[i] = 0.0;
		correct_prec[i] = 0.0;
	}

	double v = 1.0;
	for (ind_iter = pile->seqdat.begin(); ind_iter != pile->seqdat.end(); ++ind_iter) // sum over all individuals
	{
		if (ind_iter->depth == 0) // missing data for individual
			continue;

		maxlike = -1.0/0.0;

		for (j = 0; j < fnvals->optp; ++j)
		{
			for (i = 0; i < 6; ++i) // sum over all genotypic configurations
				geno_diff[i][j] = 0.0;
			indsum[j] = 0.0;
		}

		confign = 0;
		for (k1 = 0; k1 <= 2; k1 += 2) // sum over genotype1
		{
			for (k2 = 0; k2 <= 2; ++k2) // sum over genotype2
			{
				indlike[k1][k2] = 0.0;
				innerlogder[mindex] = 0.0; // d/dm P(G1,G2|f) / P(G1,G2|f)
				innerlogder[findex] = dprior[k2]; // d/df P(G1,G2|f) / P(G1,G2|f)
				v = 1.0;
				for (i = 0; i < fnvals->optp; ++i)
					dread_prob[i] = 0.0;
				for (readdat = ind_iter->rdat.begin(); readdat != ind_iter->rdat.begin() + ind_iter->depth; ++readdat) // product over all reads of individual j
				{
					epsilon = pow( 10, - (static_cast <double> (readdat->second)) / 10 );
					v = readProb(p[mindex], epsilon, k1, k2, readdat->first);
					diffRead(epsilon, k1, k2, readdat->first, &dread_prob[mindex], &dread_prob[findex]);
					indlike[k1][k2] += log(v);
					for (i = 0; i < fnvals->optp; ++i)
						innerlogder[i] += (dread_prob[i]/v) * genoprior[k2]; // P(read)'/ P(read)
				}

				if (indlike[k1][k2] > maxlike)
					maxlike = indlike[k1][k2];
				for (i = 0; i < fnvals->optp; ++i)
					geno_diff[confign][i] = innerlogder[i];
				++confign;
			}
		}
		geno_marg = 0.0;
		confign = 0;
		maxprob = exp(maxlike);
		for (k1 = 0; k1 <= 2; k1 += 2)
			for (k2 = 0; k2 <= 2; ++k2)
			{
				indl = exp(indlike[k1][k2] - maxlike);
				geno_marg += indl * genoprior[k2];
				for (i = 0; i < fnvals->optp; ++i)
					geno_diff[confign][i] *= indl * maxprob;
				++confign;
			}
		for (j = 0; j < fnvals->optp; ++j)
		{
			for (i = 0; i < 6; ++i)  // sum over all genotypic configurations
				indsum[j] += geno_diff[i][j];
			if (indsum[j] < 0.0)
				kahanSum(-exp(log(-indsum[j]) - log(geno_marg*maxprob)), &grad[j], &correct_prec[j]);
			else if (indsum[j] > 0.0)
				kahanSum(exp(log(indsum[j]) - log(geno_marg*maxprob)), &grad[j], &correct_prec[j]);
			else
				kahanSum(0.0, &grad[j], &correct_prec[j]);
		}
	}
	for (i = 0; i < fnvals->optp; ++i)
		grad[i] *= -1.0;
}

void Stats::diffRead (const double err, const int g1, const int g2, const int read, double* dm, double* df)
{
	// diffRead; derivative of the sum over true bases of P(x|b)(b|G1,G2)

	*df = 0.0;
	*dm = 0.0;
	if (g1 == 0)
	{
		if (g2 == 0)
		{
			if (read == 0)
				*dm = 0.0;
			else if (read == 1)
				*dm = 0.0;
			else
				invalidRead();
		}
		else if (g2 == 1)
		{
			if (read == 0)
				*dm = 2/3*err - 0.5;
			else if (read == 1)
				*dm = 0.5 + err/3 - err;
			else
				invalidRead();
		}
		else if (g2 == 2)
		{
			if (read == 0)
				*dm = err + err/3 - 1.0;
			else if (read == 1)
				*dm = 1.0 - err*4/3;
			else
				invalidRead();
		}
		else
			invalidG2();
	}
	else if (g1 == 2)
	{
		if (g2 == 0)
		{
			if (read == 0)
				*dm = 1.0 - err - err/3;
			else if (read == 1)
				*dm = 4/3*err - 1.0;
			else
				invalidRead();
		}
		else if (g2 == 1)
		{
			if (read == 0)
				*dm = (1.0-err)/2 - err/6;
			else if (read == 1)
				*dm = err - 0.5 - err/3;
			else
				invalidRead();
		}
		else if (g2 == 2)
		{
			if (read == 0)
				*dm = 0.0;
			else if (read == 1)
				*dm = 0.0;
			else
				invalidRead();
		}
		else
			invalidG2();
	}
	else
		invalidG1();
}

double Stats::diffGenoPrior (const double f, const int g1, const int g2)
{
	// diffGenoPrior; derivative of P(G1,G2|f)

	double dprior = 0; // d/df P(G1,G2|f)

	if (g1 == 0)
	{
		if (g2 == 0)
			dprior = -1 + f;
		else if (g2 == 1)
			dprior = 1 - 2 * f;
		else if (g2 == 2)
			dprior = f;
		else
			invalidG2();
	}
	else if (g1 == 2)
	{
		if (g2 == 0)
			dprior = -1 + f;
		else if (g2 == 1)
			dprior = 1 - 2 * f;
		else if (g2 == 2)
			dprior = f;
		else
			invalidG2();
	}
	else
		invalidG1();

	return dprior;
}
*/
