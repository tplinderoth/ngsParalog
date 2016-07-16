/*
 * stats.h
 *
 */

#ifndef STATS_H_
#define STATS_H_

#include "optim.h"

namespace Stats {
	double mafguess (Pileup* pile, bool wt);
	double calcLR (const double null, const double alt, bool islog);
	double optimLR (Optim* null, Optim* alt, double (*fn)(const double x[], const void*), void (*dfn)(const double x[], double y[], const void*), int islog, int* status);
	double negLogfn (const double para [], const void *generic_dat); // paralog likelihood function
	double readProb (const double m, const double qscore, const int g1, const int g2, const char major, const char minor, const char obs);
	void diffRead (const double err, const int g1, const int g2, const int read, double* dm, double* df);
	double genoPrior (const double f, const int g2);
	double diffGenoPrior (const double f, const int g1, const int g2);
	void kahanSum(double summand, double* total, double* comp);
	void genoErr(const char g1, const char g2);
}

#endif /* STATS_H_ */
