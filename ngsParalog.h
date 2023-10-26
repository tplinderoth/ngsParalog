//ngsParalog.h

#ifndef NGSPARALOG_H_
#define NGSPARALOG_H_

#include <fstream>
#include <iostream>
#include "Argparser.h"

const char* version = "1.3.4"; // 2023-10-26

// FUNCTION PROTOTYPES

int parseInput (int argc, char** argv, std::string& infile, std::string& outfile,
		int& quality, int& minqual, unsigned int& minind, unsigned int& mincov,
		bool& printML, int& numericGrad, int& verbose, const char* v = version);
void boundConditions (double min [], double max [], int nbounds [], int dim);
int setOptim (Optim& model, bool isalt, int verb);
void printLR (const std::string chr, const unsigned int pos, Optim* null, Optim* alt, double lr, std::ostream& os,
		const unsigned int nameBuffer, bool printml);
unsigned int calcIdBuffer (const unsigned int nameLength);
unsigned int numCovered (const Pileup* data, const unsigned int mincov);
int doLR (Argparser* params); // calculate per site LR of mismapped reads
int doRegion(); // do HMM to find coordinates of mismapping
int processPileup (std::istream& indat, std::ostream& os, Optim* altmodel, Optim* nullmodel, double (*fn)(const double x[], const void*),
		void (*dfn)(const double x[], double y[], const void*), const double qoffset, const double minq, const unsigned int mindepth,
		const unsigned int minind, int printML, int runinfo);
void endMessage (int status);

#endif /* NGSPARALOG_H_ */
