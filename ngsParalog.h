//ngsParalog.h

#ifndef NGSPARALOG_H_
#define NGSPARALOG_H_

#include <fstream>
#include <iostream>

const char* version = "0.9.1"; // version 28 July 2015

// FUNCTION PROTOTYPES
int parseInput (int argc, char** argv, std::string& infile, std::string& outfile,
		int& quality, int& minqual, unsigned int& minind, unsigned int& mincov,
		bool& printML, int& numericGrad, int& verbose, const char* v = version);
int setStreams (std::fstream& fin, std::istream& indat, std::fstream& fout, std::ostream& os, std::string* infile, std::string* outfile);
void boundConditions (double min [], double max [], int nbounds [], int dim);
bool setOptim (Optim& model, const unsigned int dim, const unsigned int npoints, int verb);
void printVal (const std::string chr, const unsigned int pos, Optim* null, Optim* alt,
		double lr, std::ostream& os, bool printml, double precision);
unsigned int nameLength (std::string name);
unsigned int numCovered (const Pileup* data, const unsigned int mincov);
void endMessage (bool error);
void info (const char*, const int qcode, const int minq, const unsigned int minind, const unsigned int mincov);
double admixOptim (FunData* optdata, Optim* alt, Optim* null, Pileup* seqdata, int numericMethod, bool logscale);
int processSeqDat (std::istream& indat, Pileup* piledat, std::ostream& os, const unsigned int mindepth, const unsigned int minInd, int numericGrad, bool printML, int verbose);
/*
int processSeqDat (std::istream& indat, Pileup* piledat, FunData* optdat, Optim* altmodel, Optim* nullmodel, std::ostream& os,
		int numericGrad, bool neglog, const int* nullindex, double increment, double decrease, int attempt, int opt_cutoff,
		double opt_thresh, Matrix<double>* altstart, Matrix<double>* nullstart, bool printML);
*/
/*
 * Arguments for processSeqDat:
 * indat - input stream who's buffer can be set to a an fstream's buffer or stdin
 * piledat - Pileup object to store sequencing data
 * optdat - FunData object that stores model information used for parameter optimization
 * altmodel - the alternative model Optim object
 * nullmodel - the null model Optim object
 * os - the output stream who's buffer can be set to an fstream's buffer or stdout
 * numericGrad - if 0 use analytic gradient, if 1 use numeric gradient in optimization
 * neglog - if 'true' likelihoods are -log scaling, if 'false' likelihoods are in probability space
 * nullindex - index in full model parameter array that are parameters used in null model
 * increment - how much to increment parameter value for each optimization attempt
 * decrease - factor to decrease 'increment' variable by in an entirely new optimization attempt
 * attempt - number of optimization tries using 'increment'/'decrease' approach to spawn starting points prior to using random start points
 * opt_cutoff - number of optimization tries before giving up moving to next site
 * opt_thresh - threshold below which likelihood ratio must be for it to be accepted
 * altstart - matrix of starting points for alternative optimization
 * nullstart - matrix of starting points for null model optimization
 * printML - if 'true' ML estimates of parameters are printed
 */

#endif /* NGSPARALOG_H_ */
