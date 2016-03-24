/*
 *  ngsParalog.cpp
 *
 *	hidden arguments (for troubleshooting):
 *  -numericGrad INT; runs the bfgs optimization using a numeric gradient if 1 [0]
 *  -verbose INT; amount of optimization output: < 0 is none, = 1 is one line, 2-98 is some, 99-100 is a lot, > 100 is max [-1]
 *  -printML 0|1; prints ML estimates of alternate allele frequency (field 6) and admixture proportion (field 7) if set to 1
*/


/*
 * TO DO:
 * 1) start null optimization with f = minor allele frequency
 * 2) start alternate optimization at null position
 * 3) fix missing data for last individual issue
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <time.h>
#include "generalUtils.h"
#include "parsePileup.h"
#include "Like.h"
#include "optim.h"
#include "stats.h"
#include "ngsParalog.h"

int main (int argc, char** argv)
{
	// declare and initialize variables
	std::string infile; // input pileup file
	std::string outfile; // output file
	int offsetqual = 33; // 33 = ASCII 33 encoding (illumina v 1.8+), 64 = ASCII 64 encoding (illumina v 1.3+)
	int minqual = 13; // minimum quality score for read to be considered
	unsigned int minind = 1; // minimum number of individuals to analyze site
	unsigned int mincov = 1; // minimum coverage for an individual to be considered "covered"
	int verbose = -1; // controls verbosity of optimization
	bool printML = false; // controls whether to print ML estimates of parameters
	int numericGrad = 0; // controls method for calculating the gradient in optimization
	bool error = false; // denotes if exiting with an error

	// print info if no arguments supplied
	if (argc==1)
	{
		info(version, offsetqual, minqual, minind, mincov);
		return 0;
	}

	// read-in, check, and assign input parameters
	int validin = parseInput (argc, argv, infile, outfile, offsetqual, minqual, minind, mincov, printML, numericGrad, verbose);
	if (validin)
		return 0; // terminate

	// set up input/output
	std::fstream fin;
	std::fstream fout;
	std::istream& indat = std::cin;
	std::ostream os(std::cout.rdbuf());
	error = setStreams(fin, indat, fout, os, &infile, &outfile);
	if (error)
	{
		endMessage(error);
		return 0;
	}

	// initialize Pileup object
	Pileup piledat;
	piledat.setQualCode(offsetqual);
	piledat.setMinQ(minqual);

	// process the data and calculate LR
	error = processSeqDat (indat, &piledat, os, mincov, minind, numericGrad, printML, verbose);
	if (error)
	{
		endMessage(error);
		return 0;
	}

	// close files
	if (infile != "-")
		fin.close();
	if (!outfile.empty())
		fout.close();

	endMessage(error);
	return 0;
}

// setStreams opens input and output stream and manages buffers
int setStreams (std::fstream& fin, std::istream& indat, std::fstream& fout, std::ostream& os, std::string* infile, std::string* outfile)
{
	// open input stream
	if (*infile != "-")
	{
		if (getFILE(fin, infile->c_str(), "in")) // open the pileup file
		{
			indat.rdbuf(fin.rdbuf()); // change input buffer from stdin to fstream's buffer if reading from file
			std::cerr << "Reading sequencing data from " << *infile << "\n";
		}
		else
			return 1;
	}
	else
		std::cerr << "Reading sequencing data from standard input\n";

	// open output stream
	if (!outfile->empty())
	{
		if (getFILE(fout, outfile->c_str(), "out"))
		{
			os.rdbuf(fout.rdbuf()); // switch output stream's buffer to the output file's buffer
			std::cerr << "Dumping results to " << *outfile << "\n";
		}
		else
			return 1;
	}
	else
		std::cerr << "Dumping results to standard output\n";
	return 0;
}

// reads-in and checks input parameters
int parseInput (int argc, char** argv, std::string& infile, std::string& outfile, int& offsetqual,
		int& minqual, unsigned int& minind, unsigned int& mincov, bool& printML, int& numericGrad, int& verbose, const char* v)
{
	int argPos = 1;
	int increment = 0;
	if ( strcmp(argv[argPos], "-help") == 0)
	{
		info(v, offsetqual, minqual, minind, mincov);
		return 1;
	}
	else
		if ((argc-1) % 2 != 0)
		{
			fprintf(stderr, "Arguments are incomplete\n--> exiting\n");
			return 1;
		}
	while (argPos < argc)
	{
		if ( strcmp(argv[argPos], "-infile") == 0 )
		{
			infile = (fexists(argv[argPos + 1]) || strcmp(argv[argPos + 1], "-") == 0) ? argv[argPos + 1] : "";
			if (infile.empty())
			{
				fprintf(stderr, "Must supply valid -infile argument\n--> exiting\n");
				return 1;
			}
		}
		else if ( strcmp(argv[argPos], "-outfile") == 0 )
			outfile = argv[argPos + 1];
		else if ( strcmp(argv[argPos], "-offsetQ") == 0 )
		{
			offsetqual = atoi(argv[argPos + 1]);
			if (offsetqual < 0)
			{
				fprintf(stderr, "-offsetQ must be a positive integer\n--> exiting\n");
				return 1;
			}
		}
		else if ( strcmp(argv[argPos], "-minQ") == 0 )
		{
			minqual = atoi(argv[argPos + 1]);
			if (minqual < 0)
			{
				fprintf(stderr, "-minQ must be a positive integer\n--> exiting\n");
				return 1;
			}
		}
		else if ( strcmp(argv[argPos], "-minind") == 0 )
		{
			minind = atoi(argv[argPos + 1]);
			if (minind < 0)
			{
				fprintf(stderr, "-minind must be a positive integer\n--> exiting\n");
				return 1;
			}
		}
		else if ( strcmp(argv[argPos], "-mincov") == 0 )
		{
			mincov = atoi(argv[argPos + 1]);
			if (mincov < 0)
			{
				fprintf(stderr, "-mincov must be a positive integer\n-->exiting\n");
				return 1;
			}
		}
		else if ( strcmp(argv[argPos], "-printML") == 0 )
		{
			printML = atoi(argv[argPos + 1]);
			if (printML != 0)
				if (printML != 1)
				{
					fprintf(stderr, "-printML must be 0 or 1\n--> exiting\n");
					return 1;
				}
		}
		else if ( strcmp(argv[argPos], "-numericGrad") == 0 )
		{
			numericGrad = atoi(argv[argPos + 1]);
			if (numericGrad != 0)
				if (numericGrad != 1)
				{
					fprintf(stderr, "-numericGrad must be 0 or 1\n--> exiting\n");
					return 1;
				}
		}
		else if ( strcmp(argv[argPos], "-verbose") == 0 )
		{
			if (argv[argPos + 1][0] == '-' && isdigit(argv[argPos + 1][1]))
				verbose = -1;
			else if ( isdigit(argv[argPos + 1][0]) )
				verbose = atoi(argv[argPos + 1]);
			else
			{
				fprintf(stderr, "Invalid argument to -verbose\n--> exiting\n");
				return 1;
			}
		}
		else
		{
			fprintf(stderr, "\nUnknown argument: %s\n\n", argv[argPos]);
			return 1; // terminate
		}

		argPos += 2 + increment;
	}

	return 0;
}

// sets values for optimization boundary conditions
void boundConditions (double min [], double max [], int nbounds [], int dim)
{
		int i = 0;

		for (i = 0; i < dim; i++)
			min[i] = 0.0;

		for (i = 0; i < dim; i++)
			max[i] = 1.0;

		for (i = 0; i < dim; i++)
			nbounds[i] = 2;
}

// sets up an optim object
bool setOptim (Optim& model,const unsigned int dim, const unsigned int npoints, int verb)
{
	// set parameter boundary conditions
	double min [dim]; // lower bounds for params
	double max [dim]; // upper bounds for params
	int nbd [dim]; // boundary constraint option for bfgs optimization
	boundConditions(min, max, nbd, dim);
	double par[dim]; // alt = {m, f}, null = {f}

	if (sizeof(par)/sizeof(par[0]) != dim) // error check
	{
		fprintf(stderr, "number of params to optimize does not equal number of provided params for alternative model\n");
		return false;
	}

	// set up optim object
	if (!model.setParN(dim))
		return false;
	if (!model.setParam(dim, par))
		return false;
	model.setVerbose(verb);
	if (!model.setBoundCntrl(max, dim, min, dim, nbd, dim))
		return false;
	model.start.allocate(npoints, dim);

	return true;
}

int processSeqDat (std::istream& indat, Pileup* piledat, std::ostream& os, const unsigned int mindepth, const unsigned int minInd, int numericGrad, bool printML, int verbose)
{

	double increment = 0.3;
	double decrease = 0.5;
	int attempt = 1;
	int opt_cutoff = 50;
	double opt_thresh = -1e-07;
	//unsigned int point_reserve = 1000;
	const unsigned int ndimensions = 2;
	bool neglog = true; // 1 = likelihoods are -log, otherwise 0

	// initialize optimization objects
	FunData optdat;
	optdat.np = ndimensions;
	const unsigned int nulldim = 1; // number params to optimize in null model
	const unsigned int naltpoints = 3; // number of optimization starting points for alternative model
	const unsigned int nnullpoints = 1; // number of optimization starting points for null model
	const double optpoints[ndimensions] = {0.2, 0.1}; // f & m anchor points for optimization
	const double optstep[ndimensions] = {0.6, 0.9}; // f & m step distance from optimization anchor
	const int nullindex[nulldim] = {0}; // index of parameters to optimize in null model
	Optim altmodel;
	if (!setOptim(altmodel, ndimensions, naltpoints, verbose))
	{
		std::cerr << "Problem encountered in function setOptim for alternate model\n--> exiting\n";
		return 0;
	}
	//altmodel.start.allocate(naltpoints, ndimensions, 0);
	Optim nullmodel;
	if (!setOptim(nullmodel, nulldim, nnullpoints, verbose))
	{
		std::cerr << "Problem encountered in function setOptim for null model\n--> exiting\n";
		return 0;
	}
	//nullmodel.start.allocate(nnullpoints, ndimensions, 0);
	Matrix<double> altstart = Optim::set2DStartVal(ndimensions, optpoints, optstep); // alt model parameter starting values for optimization
	Matrix<double> nullstart = Optim::set2DStartVal(nulldim, &optpoints[nullindex[0]], &optstep[nullindex[0]]); // null model parameter starting values for optimization

	// loop through pileup input
	std::vector<std::string>::const_iterator pline;
	std::string(line);
	Stats stat;
	stat.seed();
	double lr = 0.0;
	while (getline(indat, line))
	{
		// store data from pileup line in Pileup object
		piledat->getSeqDat(line);
		if (piledat->fail())
		{
			fprintf(stderr, "Error in parsing sequencing data - check pileup\n");
			return 1;
		}
		if (numCovered(piledat, mindepth) < minInd)
		{
			fprintf(stderr, "skipping %s %u --> inadequate coverage\n", piledat->seqName().c_str(), piledat->position());
			continue;
		}
		optdat.sitedat = piledat;
		// perform optimization and calculate LR
		lr = admixOptim(&optdat, &altmodel, &nullmodel, piledat, numericGrad, neglog);
		if (nullmodel.fail() || altmodel.fail())
			return 1;
		if (lr < opt_thresh)
			lr = stat.optimLR(Like::negLogfn, Like::calcGradient, &optdat, &altmodel, &nullmodel, numericGrad, neglog, nullindex, increment, decrease, attempt, opt_cutoff, opt_thresh, &altstart, &nullstart, verbose);
		if (stat.fail())
		{
			fprintf(stderr, "Error occurred within Stats::optimLR\n");
			return 1;
		}
		// dump results
		printVal(piledat->seqName(), piledat->position(), &nullmodel, &altmodel, lr, os, printML, opt_thresh);
	}
	return 0;
}

double admixOptim (FunData* optdata, Optim* alt, Optim* null, Pileup* seqdata, int numericMethod, bool logscale)
{
	// optimize null
	optdata->optp = null->getDim();
	null->start[0][0] = seqdata->altfreq();
	null->runMultiOptim(Like::negLogfn, Like::calcGradient, optdata, numericMethod);
	if (null->fail())
	{
		fprintf(stderr, "Error in optimization for %s position %u\n", seqdata->seqName().c_str(), seqdata->position());
		return -1.0/0.0;
	}
	// optimize alternative
	optdata->optp = alt->getDim();
	alt->start[0][0] = null->getMLParam(0);
	alt->start[0][1] = 1.0;
	static double mvals [] = {0.5, 0.1};
	static double mlen = sizeof(mvals)/sizeof(double);
	double f = 0.0;
	for (int i = 0; i < mlen; ++i)
	{
		alt->start[i+1][1] = mvals[i];
		f = (1/alt->start[i+1][1]) * seqdata->altfreq();
		alt->start[i+1][0] = f < 1.0 ? f : 1.0;
	}
	alt->runMultiOptim(Like::negLogfn, Like::calcGradient, optdata, numericMethod);
	if (alt->fail())
	{
		fprintf(stderr, "Error in optimization for %s position %u\n", seqdata->seqName().c_str(), seqdata->position());
		return -1.0/0.0;
	}
	// calculate LR
	return (Stats::calcLR(null->getNegllh(), alt->getNegllh(), logscale));
}

unsigned int numCovered (const Pileup* data, const unsigned int mincov)
{
	static std::vector<SiteData>::const_iterator indIter;
	unsigned int ncovered = 0;
	for (indIter = data->seqdat.begin(); indIter != data->seqdat.end(); ++indIter)
	{
		if (indIter->depth >= mincov)
			++ncovered;
	}
	return ncovered;
}

// printVal prints output
void printVal (const std::string chr, const unsigned int pos, Optim* null, Optim* alt, double lr, std::ostream& os, bool printml, double precision)
{
	static unsigned int namesz = 0;
	static int i = 0;
	if (namesz == 0)
		namesz = nameLength(chr);
	if (lr < 0.0) // prevents reporting negative zero - for output aesthetics
	{
		if (lr > precision)
			lr *= -1.0;
	}
	os << std::setw(namesz) << std::left << chr << '\t' << std::setw(12) << pos << '\t' << std::fixed << std::setw(12) << std::right
			<< null->getNegllh() << '\t' << std::setw(12) << std::right << alt->getNegllh() << '\t' << std::setw(12) << std::right << lr;
	if (printml)
		for (i = 0; i < alt->getDim(); ++i)
			os << '\t' << std::setw(12) << std::right << alt->getMLParam(i);
	os << '\n';
}

unsigned int nameLength (std::string name)
{
	int extra = 9;
	unsigned int len = name.length() + extra;
	while (len%8)
		++len;
	return len;
}

// info prints information about the program
void info (const char* v, const int qcode, const int minq, const unsigned int minind, const unsigned int mincov)
{
	int w1 = 10;
	fprintf(stderr, "\nngsParalog version %s\n", v);
	std::cerr << "\nInput:"
	<< "\n" << std::setw(w1) << std::left << "-infile" << std::setw(w1) << "FILE|-" << "pileup format file of reads and quality scores; specifying '-' will read from STDIN [NULL]"
	<< "\n" << std::setw(w1) << std::left << "-outfile" << std::setw(w1) << "FILE" <<  "name of output file; if not provided results printed to STDOUT [NULL]"
	<< "\n" << std::setw(w1) << std::left << "-minQ" << std::setw(w1) << "INT" << "minimum quality score for a read to be kept [" << minq << "]"
	<< "\n" << std::setw(w1) << std::left << "-mincov" << std::setw(w1) << "INT" << "minimum per individual coverage for an individual to be considered covered for -minind [" << mincov << "]"
	<< "\n" << std::setw(w1) << std::left << "-minind" << std::setw(w1) << "INT" << "minimum number of individuals with -mincov coverage for a site to be analyzed [" << minind << "]"
	<< "\n" << std::setw(w1) << std::left << "-offsetQ" << std::setw(w1) << "INT" << "minimum possible ASCII decimal value used to encode quality scores [" << qcode << "]"
	<< "\n\nOutput by field:"
	<< "\n(1) chromosome name"
	<< "\n(2) site"
	<< "\n(3) -log null likelihood"
	<< "\n(4) -log alternative likelihood"
	<< "\n(5) likelihood ratio"
	<< "\n\n";
}

void endMessage (bool error)
{
	if (!error)
		std::cerr << "Finished!\n";
	else
		std::cerr << "-->exiting\n";
}
