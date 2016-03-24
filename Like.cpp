// Like.cpp

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include "optim.h"
#include "Like.h"
#include "parsePileup.h"
//#include <iostream> // debug

// negLogfn is the main likelihood function
double Like::negLogfn (const double para [], const void *generic_dat)
{
	// para[0] = alternate allele frequency; para[1] = admixture proportion

	const FunData* fnvals = static_cast<const FunData*>(generic_dat);
	const Pileup* pile = static_cast<const Pileup*>(fnvals->sitedat);
	const int npara = 2; // number parameters in model
	const int findex = 0; // index position of allele frequency parameter
	const int mindex = npara - 1 - findex; // index position of admixture proportion
	static double p [npara];
	static double epsilon; // Phred scaled quality score
	static double geno_marg; // marginal sum over genotype2 configurations
	static double negloglike; // negative log likelihood of likelihood function
	static double genoprior [3];
	static int k1, k2;
	static Matrix<double> indlike (3, 3); // llh for individual at all genotypic configurations
	static double maxlike; // max llh for individual among all genotypic configurations
	static int isnull;
	static double compensation = 0.0; // compensation for lost low-order bits in Kahan Sum

	isnull = (fnvals->optp < fnvals->np) ? 1 : 0;
	if (isnull == 1) // null case; fix m to 1
	{
		p[findex] = para[0]; // alternate allele frequency (f)
		p[mindex] = 1.0; // admixture proportion (m)
	}
	else // alternative case
	{
		for (int i = 0; i < fnvals->np; ++i)
				p[i] = para[i];
	}
	compensation = 0.0;
	negloglike = 0.0;

	// P(G1) does not depend on f; k1 can be 0 or 2
	genoprior[0] = genoPrior(p[findex], 0, 0); // P(G1 = {0,2}, G2 = 0|f)
	genoprior[1] = genoPrior(p[findex], 0, 1);  // P(G1 = {0,2}, G2 = 1|f)
	genoprior[2] = genoPrior(p[findex], 0, 2); // P(G1 = {0,2}, G2 = 2|f)
	static std::vector<SiteData>::const_iterator ind_iter;
	static std::vector< std::pair<int, int> >::const_iterator readdat;

	for (ind_iter = pile->seqdat.begin(); ind_iter != pile->seqdat.end(); ++ind_iter) // sum over all individuals
	{
		if (ind_iter->depth == 0) // missing data for individual
			continue;
		maxlike = -1.0/0.0;
		for (k1 = 0; k1 <= 2; k1 += 2) // sum over genotype1
		{
			for (k2 = 0; k2 <= 2; ++k2) // sum over genotype2
			{
				indlike[k1][k2] = 0.0;
				for (readdat = ind_iter->rdat.begin(); readdat != ind_iter->rdat.begin() + ind_iter->depth; ++readdat) // product over all reads of individual
				{
					epsilon = pow( 10, - (static_cast <double > (readdat->second)) / 10 );
					indlike[k1][k2] += log( readProb(p[mindex], epsilon, k1, k2, readdat->first) );
				}
				if (indlike[k1][k2] > maxlike)
					maxlike = indlike[k1][k2];
			}
		}
		geno_marg = 0.0;
		for (k1 = 0; k1 <= 2; k1 += 2)
			for (k2 = 0; k2 <= 2; ++k2)
				geno_marg += exp(indlike[k1][k2] - maxlike) * genoprior[k2]; // switch back to probability space
		kahanSum(-(log(geno_marg) + maxlike), &negloglike, &compensation);
	}
	return negloglike;
}

// calcGradient gets the gradient of the negative log likelihood function
void Like::calcGradient (const double para [], double grad [], const void* generic_dat)
{
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

// readProb; sum over true bases of P(x|b)(b|G1,G2)
double Like::readProb (const double m, const double err, const int g1, const int g2, const int read)
{
	double pread = 1; // probability of read
	if (g1 == 0)
	{
		if (g2 == 0)
		{
			if (read == 0)
				pread = 1.0 - err;
			else if (read == 1)
				pread = err;
			else
				invalidRead();
		}
		else if (g2 == 1)
		{
			if (read == 0)
				pread = 1.0 - err + m*(2*err/3 - 0.5);
			else if (read == 1)
				pread = err + m*(0.5 + err/3) - err*m;
			else
				invalidRead();
		}
		else if (g2 == 2)
		{
			if (read == 0)
				pread = 1.0 - err + m*(err + err/3 - 1.0);
			else if (read == 1)
				pread = err + m*(1.0 + (2*err)/3 - 2*err);
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
				pread = m*(1.0 - err - err/3) + err/3;
			else if (read == 1)
				pread = 2*err*(m + 1/3 - m/3) - m + 1.0 - err;
			else
				invalidRead();
		}
		else if (g2 == 1)
		{
			if (read == 0)
				pread = 2*m*(0.25 - err/3) + err/3;
			else if (read == 1)
				pread = err*(m + 2/3 - m/3 - 1.0) + m/2 + 1.0 - m;
			else
				invalidRead();
		}
		else if (g2 == 2)
		{
			if (read == 0)
				pread = err/3;
			else if (read == 1)
				pread = 1.0 - err/3;
			else
				invalidRead();
		}
		else
			invalidG2();
	}
	else
		invalidG1();

	return pread;
}

// diffRead; derivative of the sum over true bases of P(x|b)(b|G1,G2)
void Like::diffRead (const double err, const int g1, const int g2, const int read, double* dm, double* df)
{
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

// genoPrior; P(G1,G2|f)
double Like::genoPrior (const double f, const int g1, const int g2)
{
	double prob = 0;

	if (g1 == 0)
	{
		if (g2 == 0)
			prob = 0.5 * (1 - f) * (1 - f);
		else if (g2 == 1)
			prob = f * (1 - f);
		else if (g2 == 2)
			prob = 0.5 * f * f;
		else
			invalidG2();
	}
	else if (g1 == 2)
	{
		if (g2 == 0)
			prob = 0.5 * (1 - f) * (1 - f);
		else if (g2 == 1)
			prob = f * (1 - f);
		else if (g2 == 2)
			prob = 0.5 * f * f;
		else
			invalidG2();
	}
	else
		invalidG1();

	return prob;
}

// diffGenoPrior; derivative of P(G1,G2|f)
double Like::diffGenoPrior (const double f, const int g1, const int g2)
{
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

void Like::kahanSum(double summand, double* total, double* comp)
{
	double x = summand - *comp;
	double y = *total + x;
	*comp = (y - *total) - x;
	*total = y;
}

void Like::invalidRead ()
{
	fprintf(stderr, "\nInvalid read type encountered\n--> exiting\n");
	exit(EXIT_FAILURE);
}

void Like::invalidG1 ()
{
	fprintf(stderr, "\nInvalid number of non-ref alleles for locus 1 encountered\n--> exiting\n");
	exit(EXIT_FAILURE);
}

void Like::invalidG2 ()
{
	fprintf(stderr, "\nInvalid number of non-ref alleles for locus 2 encountered\n--> exiting\n");
	exit(EXIT_FAILURE);
}
