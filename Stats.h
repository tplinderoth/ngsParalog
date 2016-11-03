/*
 * stats.h
 *
 */

#ifndef STATS_H_
#define STATS_H_

#include "optim.h"

class readProb {
	typedef double (*ReadPrFn) (double m, double err);
public:
	// public member functions
	readProb();
	// public member variables
	ReadPrFn majorprobs [3]; // matrix of major read probability functors
	ReadPrFn minorprobs [3]; // matrix of minor read probability functors
};

namespace Stats {
double mafguess (Pileup* pile, bool wt);
double calcLR (const double null, const double alt, bool islog);
double optimLR (Optim* null, Optim* alt, double (*fn)(const double x[], const void*), void (*dfn)(const double x[], double y[], const void*), int islog, int* status);
double negLogfn (const double para [], const void *generic_dat); // paralog likelihood function
double prRead (double m, double qscore, unsigned int g2, char major, char minor, char obs);
void diffRead (const double err, const int g1, const int g2, const int read, double* dm, double* df);
double genoPrior (const double f, const int g2);
double diffGenoPrior (const double f, const int g1, const int g2);
void kahanSum(double summand, double* total, double* comp);
void genoErr(unsigned int g1, unsigned int g2);
// major read probabilities
double major00 (double m, double err);
double major01 (double m, double err);
double major02 (double m, double err);
//double major20 (double m, double err);
//double major21 (double m, double err);
//double major22 (double m, double err);
// minor read probabilities
double minor00 (double m, double err);
double minor01 (double m, double err);
double minor02 (double m, double err);
//double minor20 (double m, double err);
//double minor21 (double m, double err);
//double minor22 (double m, double err);
};

#endif /* STATS_H_ */
