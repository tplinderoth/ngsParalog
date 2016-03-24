/*
 * stats.cpp
 */

#include <cstdlib>
#include <math.h>
#include "optim.h"
#include "stats.h"
#include <time.h>
#include "parsePileup.h"
//#include <iostream> // debug

// stats class constructor
Stats::Stats()
	: _fail(false),
	  _step_full(0),
	  _step_reduced(0),
	  _start_full(0),
	  _start_reduced(0),
	  _step_full_size(0),
	  _step_reduced_size(0),
	  _start_full_size(0),
	  _start_reduced_size(0),
	  _stat(0),
	  _seed(0)
{}

Stats::~Stats()
{
	if (_step_full)
	{
		delete [] _step_full;
		_step_full = 0;
		_step_full_size = 0;
	}
	if (_step_reduced)
	{
		delete [] _step_reduced;
		_step_reduced = 0;
		_step_reduced_size = 0;
	}
	if (_start_full)
	{
		delete [] _start_full;
		_start_full = 0;
		_start_full_size = 0;
	}
	if (_start_reduced)
	{
		delete [] _start_reduced;
		_start_reduced = 0;
		_start_reduced_size = 0;
	}
}

void Stats::setFullStepSize (int n)
{
	if (_step_full)
		delete [] _step_full;
	_step_full = new double [n];
	_step_full_size = n;
}

void Stats::setReducedStepSize (int n)
{
	if (_step_reduced)
		delete [] _step_reduced;
	_step_reduced = new double[n];
	_step_reduced_size = n;
}

void Stats::setFullStartSize (int n)
{
	if (_start_full)
		delete [] _start_full;
	_start_full = new double [n];
	_start_full_size = n;
}

void Stats::setReducedStartSize (int n)
{
	if (_start_reduced)
		delete [] _start_reduced;
	_start_reduced = new double [n];
	_start_reduced_size = n;
}

size_t Stats::getFullStepSize ()
{
	return _step_full_size;
}

size_t Stats::getReducedStepSize ()
{
	return _step_reduced_size;
}

size_t Stats::getFullStartSize ()
{
	return _start_full_size;
}

size_t Stats::getReducedStartSize ()
{
	return _start_reduced_size;
}

// stats::optimLR calculates a likelihood ratio using MLEs for the alt and null models of <= 2 parameters
double Stats::optimLR ( double (*fun)(const double x[], const void*), void (*dfun)(const double x[], double y[], const void*),
		FunData* likedat, Optim* altopt, Optim* nullopt, int gradMeth, bool logscale, const int* reduced_ind, double increment,
		double decrease, int attempt, int tries, double lr_threshold, Matrix<double>* alt_start, Matrix<double>* null_start, int verb)
{
	static double likeratio = 0.0;
	if (!_step_full)
	{
		setFullStepSize(altopt->dim);
		setReducedStartSize(nullopt->dim);
		setReducedStepSize(nullopt->dim);
	}
	if (!_seed)
		seed();

	likeratio = multiLR(fun, dfun, likedat, altopt, nullopt, gradMeth, logscale, alt_start, null_start);
	if (likeratio < lr_threshold)
	{
		if (gradMeth != 1)
			likeratio = multiLR(fun, dfun, likedat, altopt, nullopt, 1, logscale, alt_start, null_start);

		int i = 0;
		int t = 0;
		while (!_fail && likeratio < lr_threshold)
		{
			// spawn new optimization starting points
			if (verb > 0)
				fprintf(stderr, "\nTrying new random optimization start points\n");
			for (i = 0; i < altopt->dim; ++i)
			{
				if (t <= attempt)
					_step_full[i] = increment;
				else
					_step_full[i] = altopt->lowb[i] + static_cast <double> (rand())/(RAND_MAX/(altopt->upb[i] - altopt->lowb[i]));
			}
			for (i = 0; i < nullopt->dim; ++i)
			{
				_start_reduced[i] = altopt->lowb[reduced_ind[i]];
				_step_reduced[i] = _step_full[reduced_ind[i]];
			}
			++t;
			increment *= decrease;
			Matrix<double> newalt = Optim::set2DStartVal(altopt->dim, &(*altopt->lowb), &(*_step_full));
			Matrix<double> newnull = Optim::set2DStartVal(nullopt->dim, &(*_start_reduced), &(*_step_reduced)); // alt model parameter starting values for optimization
			/*
			if ( !Optim::mod2DStartVal(altopt->dim, &(*altopt->lowb), &(*_step_full), &altopt->start) )
				_fail = true;
			if (!Optim::mod2DStartVal(nullopt->dim, &(*_start_reduced), &(*_step_reduced), &nullopt->start))
				_fail = true;
			*/
			// reoptimize
			likeratio = multiLR(fun, dfun, likedat, altopt, nullopt, gradMeth, logscale, &newalt, &newnull);
			if (gradMeth != 1 && likeratio < lr_threshold)
				likeratio = multiLR(fun, dfun, likedat, altopt, nullopt, 1, logscale, &newalt, &newnull);
			if (t >= tries)
			{
				optFailMsg(likedat, tries);
				break;
			}
		}
	}
	return likeratio;
}

// stats::calcLR calcLR calculates the likelihood ratio
double Stats::calcLR (const double null, const double alt, bool islog)
{
	double lr = -1.0/0.0;

	if (islog) // likelihoods are -log scaled
		lr = 2 * null - 2 * alt;
	else
		lr = -2 * log(null) + 2 * log(alt);
	return lr;
}

double Stats::multiLR (double (*fun)(const double x[], const void*), void (*dfun)(const double x[], double y[], const void*),
		FunData* likedat, Optim* altopt, Optim* nullopt, int gradMeth, bool logscale,
		Matrix<double>* alt_start, Matrix<double>* null_start)
{
	likedat->optp = altopt->dim;
	altopt->runMultiOptim(fun, dfun, likedat, gradMeth, alt_start);
	likedat->optp = nullopt->dim;
	nullopt->runMultiOptim(fun, dfun, likedat, gradMeth, null_start);
	return (calcLR(nullopt->negllh, altopt->negllh, logscale));
}

void Stats::optFailMsg (FunData* data, const int thresh)
{
	Pileup* info = static_cast<Pileup*>(data->sitedat);
	fprintf(stderr, "%i optimization failures for %s %u --> moving on...\n",
			thresh, (info->seqName()).c_str(), info->position());
}

// returns _fail member
bool Stats::fail ()
{
	return _fail;
}

void Stats::seed()
{
	time_t t = time(NULL);
	srand((unsigned)t);
	_seed = t;
}

time_t Stats::getSeed()
{
	return _seed;
}
