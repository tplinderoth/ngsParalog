/*
 *  ngsParalog.cpp
 *
 * ngsParalog - Genomic duplication detection from NGS data
 * Copyright (C) 2016 Tyler P. Linderoth
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Tyler's contact information: tylerp.linderoth@gmail.com
 *
 *
 *	hidden arguments (for troubleshooting):
 *  -numericGrad INT; runs the bfgs optimization using a numeric gradient if 1 [0]
 *  -verbose INT; amount of optimization output: < 0 is none, = 1 is one line, 2-98 is some, 99-100 is a lot, > 100 is max [-1]
 *  -printML 0|1; prints ML estimates of alternate allele frequency (field 6) and admixture proportion (field 7) if set to 1
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <time.h>
#include "generalUtils.h"
#include "parsePileup.h"
#include "optim.h"
#include "Stats.h"
#include "ngsParalog.h"

int main (int argc, char** argv)
{
	int rc = 0; // store function return values
	Argparser runpar; // stores input parameters

	// parse user input
	if((rc=runpar.parseInput(argc, argv, version)) != 0)
	{
		if (rc < 0)
		{
			endMessage(rc);
			return 1;
		}
		else
			return 0;
	}
	else
	{
		// set IO streams
		if((rc=runpar.setStreams(runpar.infile_name(), runpar.outfile_name())))
		{
			endMessage(rc);
			return 1;
		}
	}

	// run the task
	if (strcmp(runpar.task(), "calcLR") == 0)
		rc = doLR(&runpar);
	else if (strcmp(runpar.task(), "findRegion") == 0)
		rc = doRegion ();
	else
	{
		std::cerr << "Undefined task in call to main\n";
		rc=1;
		endMessage(rc);
	}

	// close IO streams
	runpar.closeStreams();

	// exit program
	endMessage(rc);
	return rc;
}

int doLR (Argparser* params)
{
	int rv = 0;

	// set likelihood functions
	//double (*likefn)(const double x[], const void*) = Stats::negLogfn;
	double (*likefn)(const double x[], const void*) = &Stats::negLogfn;
	void (*dlikefn)(const double x[], double y[], const void*) = NULL; // could adjust this with Argparser _numericGrade member

	// initialize optimization objects
	Optim altmodel;
	if (setOptim(altmodel, 1, params->verblevel()))
	{
		std::cerr << "Problem encountered in function setOptim for alternate model\n";
		return 1;
	}

	Optim nullmodel;
	if (setOptim(nullmodel, 0, params->verblevel()))
	{
		std::cerr << "Problem encountered in function setOptim for null model\n";
		return 1;
	}

	// analyze pileup data
	try
	{
		rv = processPileup(params->input(), params->output(), &altmodel, &nullmodel, likefn, dlikefn, params->offsetQ(),
				params->minQ(), params->mincov(), params->minind(), params->printML(), params->ploidy());
	}
	catch (PileupException& error)
	{
		std::cerr << error.what() << "\nGracefully terminating\n";
		return -1;
	}
	catch (std::exception& error)
	{
		std::cerr << error.what() << "\nGracefully terminating\n";
		return -1;
	}

	return rv;
}

int doRegion ()
{
	std::cerr << "Dang, findRegion still not implemented ...\n";
	return 0;
}

int setOptim (Optim& model, bool isalt, int verb)
{
// sets up optim objects

	const unsigned int dim = isalt ? 2 : 1; // number of parameters to optimize
	int optconditions = 6; // number of conditions that must be set to perform bfgs optimization
	const double startpoints[2] = {0.2, 0.05}; // f & m anchor points for optimization
	const double step[2] = {0.7, 0.48}; // f & m step distance from optimization anchor
	int nullidx = 0; // index of null parameters in parameter arrays
	double* min = new double[dim]; // minimum parameter values
	double* max = new double[dim]; // maximum parameter values
	unsigned int* nbd = new unsigned int[dim]; // number of boundary conditions
	unsigned int i = 0;

	// set boundary values
	for (i = 0; i < dim; i++)
		min[i] = 0.0;
	for (i = 0; i < dim; i++)
		max[i] = 1.0;
	for (i = 0; i < dim; i++)
		nbd[i] = 2;

	// set values for optim
	if (!model.setParN(dim))
		return false;
	if (!model.setBoundCntrl(max, min, nbd, dim))
		return false;

	// set amount of output
	model.setVerbose(verb);

	// set optimization start points
	if (isalt)
	{
		if (model.initStartMatrix(dim, startpoints, step, 1))
			return false;
	}
	else
	{
		if (model.initStartMatrix(dim, &startpoints[nullidx], &step[nullidx]))
			return false;
	}

	// make sure everything has been initialized
	if (model.conditions() < optconditions)
	{
		fprintf(stderr,"Optim object has only %i out of %i optimization conditions set in call to setOptim\n", model.conditions(), optconditions);
		return false;
	}

	delete [] min;
	delete [] max;
	delete [] nbd;

	return 0;
}

int processPileup (std::istream& indat, std::ostream& os, Optim* altmodel, Optim* nullmodel, double (*fn)(const double x[], const void*),
		void (*dfn)(const double x[], double y[], const void*), const double qoffset, const double minq, const unsigned int mindepth,
		const unsigned int minind, int printML, const int ploidy)
{
	const char delimiter = '\t'; // assume pileup is tab delimited
	const bool weightcount = true;
	const int neglog = 1; // 1 = likelihoods are -log, otherwise 0
	double lr = 0.0;
	int optfail = 0;
	std::string(line);

	// initialize Pileup object
	Pileup piledat;
	piledat.setQualCode(qoffset);
	piledat.setMinQ(minq);
	getline(indat, line);
	piledat.setn(line, delimiter);
  piledat.setPloidy(ploidy);

	// calculate sequence ID buffer size for pretty printing
	unsigned int idbuffer = calcIdBuffer(Pileup::idSize(line, delimiter));

	// check that there at at least minind individuals in pileup
	if (minind > piledat.nInd())
	{
		fprintf(stderr,"Fewer than %u individuals in dataset --> decrease -minind argument\n",minind);
		return 1;
	}

	// loop through pileup input
	while (!line.empty())
	{
		// store data in Pileup object
		piledat.getSeqDat(line, delimiter);

		if (piledat.fail())
		{
			fprintf(stderr, "Error parsing pileup data\n");
			return 1;
		}

		// check for excessive missing data
		if (numCovered(&piledat, mindepth) < minind)
		{
			fprintf(stderr, "skipping %s %u --> inadequate coverage\n", piledat.seqName().c_str(), piledat.position());
			getline(indat,line);
			continue;
		}

		// set major and minor allele and store pileup data in optim objects
		piledat.setMajor(piledat.empiricalMajor(weightcount)); // assign most common allele as major
		piledat.setMinor(piledat.empiricalMinorFast(weightcount)); // assign second most common allele as minor
		altmodel->data = nullmodel->data = &piledat;

		// perform optimization and calculate LR
		try
		{
			lr = Stats::optimLR(nullmodel, altmodel, fn, dfn, neglog, &optfail);
			printLR(piledat.seqName(), piledat.position(), nullmodel, altmodel, lr, os, idbuffer, printML);
		}
		catch (const OptimFailureException& error)
		{
			std::cerr << error.what() << " -> skipping " << piledat.seqName() << " " << piledat.position() << "\n";
		}
		catch (const NoDataException& error)
		{
			std::cerr << error.what() << " -> skipping site\n";
		}
		// fetch next line
		getline(indat,line);
	}

	return 0;
}

unsigned int numCovered (const Pileup* data, const unsigned int mincov)
{
	static std::vector<SiteData>::const_iterator indIter;
	unsigned int ncovered = 0;
	for (indIter = data->seqdat.begin(); indIter != data->seqdat.end(); ++indIter)
	{
		if (indIter->cov() >= mincov)
			++ncovered;
	}
	return ncovered;
}


void printLR (const std::string chr, const unsigned int pos, Optim* null, Optim* alt, double lr, std::ostream& os,
		const unsigned int nameBuffer, bool printml)
{
	const double thresh = -1e-07;
	const int prec = 8;
	int i = 0;
	if (lr < 0.0) // prevents reporting negative zero (numerical noise)
	{
		if (lr > thresh)
			lr = 0.0;
	}
	os << std::setw(nameBuffer) << std::left << chr
			<< '\t' << std::setw(12) << pos
			<< '\t' << std::fixed << std::setprecision(prec) << null->llh()
			<< '\t' << std::fixed << std::setprecision(prec) << alt->llh()
			<< '\t' << std::fixed << std::setprecision(prec) << lr;
	if (printml)
		for (i = 0; i < alt->getDim(); ++i)
			os << '\t' << std::fixed << std::setprecision(prec) << alt->mlparam(i);
	os << '\n';
}

unsigned int calcIdBuffer (const unsigned int nameLength)
{
	const int extra = 9;
	unsigned int len = nameLength + extra;
	while (len%8)
		++len;
	return len;
}

void endMessage (int status)
{
	if (!status)
		std::cerr << "Finished!\n";
	else
		std::cerr << "-->exiting\n";
}
