// optim.cpp
// last edit: 2/17/2015

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <time.h>
#include <math.h>
#include "optim.h"
#include "bfgs.h"
// #include <iostream> // debug

// Optim class constructor
Optim::Optim (int noise)
	: par(0),
	  mlpar(0),
	  upb(0),
	  lowb(0),
	  dim(0),
	  nbounds(0),
	  xtradat(0),
	  _fail(false),
	  _nconditions(0),
	  negllh( -1.0 ),
	  _failcount(0),
	  _randset(0)
{
	setVerbose(noise);
}

// Optim class destructor
Optim::~Optim()
{
	// clear parameters
	if (par)
	{
		delete [] par;
		--_nconditions;
	}
	if (mlpar)
		delete [] mlpar;

	//clear bounds
	if (nbounds)
	{
		delete [] nbounds;
		--_nconditions;
	}
	if (upb)
	{
		delete [] upb;
		--_nconditions;
	}
	if (lowb)
	{
		delete [] lowb;
		--_nconditions;
	}
}

// Optim::cpyConditions copies existing optimization info when dimension is changed
template< class T > void Optim::cpyConditions (int olddim, int newdim, T* oldarr)
{
	int arrsize;
	if (newdim < olddim)
		arrsize = newdim * sizeof(T);
	else
		arrsize = olddim * sizeof(T);
	T newarr[newdim];
	memcpy(newarr, oldarr, arrsize);
	delete [] oldarr;
	oldarr = new T[newdim];
	memcpy(oldarr, newarr, arrsize);
	if (newdim > olddim)
	{
		for (int i = olddim; i < newdim; i++)
			oldarr[i] = 0;
	}
}

// getParN gets number of parameters to optimize
bool Optim::setParN (int n)
{
	if (n > 0)
	{
		if (dim)
		{
			if (par)
				cpyConditions(dim, n, par);
			if (lowb)
				cpyConditions(dim, n, lowb);
			if (upb)
				cpyConditions(dim, n, upb);
			if (nbounds)
				cpyConditions(dim, n, nbounds);
			if (mlpar)
				cpyConditions(dim, n, mlpar);
		}
		dim = n;
	}
	else
	{
		fprintf(stderr, "Invalid number of parameters in Optim::setParN\n");
		dim = 0;
		return false;
	}
	return true;
}

// Optim::setBoundCntrl sets boundary conditions for optimization
bool Optim::setBoundCntrl (const double upperb [], int usize, const double lowerb [], int lsize, const int constraint [], int bsize)
{
	if (!setBoundCon(constraint, bsize))
		return false;
	if (!setBound(upperb, usize, lowerb, lsize))
			return false;
	return true;
}

// Optim::setUpBound sets the upper bound on parameters to optimize
bool Optim::setUpBound (const double u [], const int usize)
{
	if (usize == dim)
	{
		if (!upb)
		{
			upb = new double [dim];
			++_nconditions;
		}
		for (int i = 0; i < dim; ++i)
			upb[i] = u[i];
	}
	else
	{
		fprintf(stderr, "Number of bounds doesn't match Optim::dim in Optim::setUpBound\n");
		upb = NULL;
		return false;
	}
	return true;
}

// Optim::setLowBound sets the lower bound on parameters to optimize
bool Optim::setLowBound (const double lo [], const int lsize)
{
	if (lsize == dim)
	{
		if (!lowb)
		{
			lowb = new double [dim];
			++_nconditions;
		}
		for (int i = 0; i < dim; ++i)
			lowb[i] = lo[i];

	}
	else
	{
		fprintf(stderr, "Number of bounds doesn't match Optim::dim in Optim::setLowBound\n");
		lowb = NULL;
		return false;
	}
	return true;
}

void Optim::setMLPar()
{
	if (dim > 0)
		if (!mlpar)
        {
			mlpar = new double[dim];
            for (int i = 0; i < dim; ++i)
            	mlpar[i] = 0.0;
        }
        else
        	fprintf(stderr, "Optim::mlpar already set in call to Optim::setMLPar()\n");
	else
		fprintf(stderr, "warning: Optim::dim not set in call to Optim::setMLPar\n");
}


// Optim::setUpBound sets the upper and lower boundaries on parameters to optimize
bool Optim::setBound (const double upb [], int usize, const double lob [], int lsize)
{
	if (!setUpBound(upb, usize)) // set upper bound
		return false;
	if (!setLowBound(lob, lsize)) // set lower bound
		return false;
	return true;
}

// Optim::setBoundCon sets the number of boundaries for each parameter
bool Optim::setBoundCon (const int btype [], const int bsize)
{
	if (bsize == dim)
	{
		if (!nbounds)
		{
			nbounds = new int [dim];
			++_nconditions;
		}
		for (int i = 0; i < dim; ++i)
			nbounds[i] = btype[i];
	}
	else
	{
		fprintf(stderr, "Number of bound constraint types doesn't match number of parameters in Optim::setBoundCon\n");
		nbounds = NULL;
		return false;
	}
	return true;
}

// Optim::setVerbose sets the amount of output during optimization
void Optim::setVerbose (int v)
{
	verb = v;
}

// Optim::setExtraDat sets data used by likelihood function
void Optim::setExtraDat (void* data)
{
	if (data)
		xtradat = data;
	else
	{
		fprintf(stderr, "Assignment of NULL pointer in Optim::setExtraDat\n");
		xtradat = NULL;
	}

}

// Optim::setParam sets user defined starting parameter values
bool Optim::setParam (const int npars, double* vals)
{
	if (npars == dim)
	{
		if (!par)
		{
			par = new double [dim];
			++_nconditions;
		}
		for (int i = 0; i < dim; ++i)
			par[i] = vals[i];
		if (!mlpar)
			setMLPar();
	}
	else
	{
		fprintf(stderr, "Number of parameters supplied to Optim::setParam differs from Optim::dim\n");
		par = NULL;
		return false;
	}
	return true;
}

// Optim::setRandParam sets random starting parameter values in [0,1]
void Optim::setRandParam (const int npars)
{
	if (!_randset)
		seed();
	if (npars == dim)
	{
		if (!par)
		{
			par = new double [dim];
			++_nconditions;
		}
		for (int i = 0; i < dim; ++i)
			par[i] = (static_cast <double> (rand()) / (RAND_MAX) );
	}
	else
	{
		fprintf(stderr, "Number of parameters supplied to Optim::setRandParam differs from Optim::dim\n");
		par = NULL;
	}
}

void Optim::seed()
{
	time_t t = time(NULL);
	srand((unsigned)t);
	_randset = t;
}

time_t Optim::getSeed()
{
	return _randset;
}

// Optim::runBFGS calls the L-BFGS-B optimization
double Optim::runBFGSB (double (*fn)(const double x[], const void*), void (*dfn)(const double x[], double y[], const void*))
{
	/*
	run bounded BFGS
	fn is function to be minimized ( -log() )
	dfn is the function to the calculate the gradient of the function
	can use 'NULL' for gradient function to calculate numerical gradient
	*/
	if (_nconditions < 4)
	{
		fprintf(stderr, "All necessary optimization conditions not set in Optim::runBFGS\n");
		_fail = true;
		return negllh = -1.0;
	}
	negllh = findmax_bfgs(dim, par, xtradat, fn, dfn, lowb, upb, nbounds, verb, &_failcount);
	return negllh;
}

// Optim::clearBounds frees memory allocated for bounding
void Optim::clearBounds ()
{
	if (nbounds)
	{
		delete [] nbounds;
		nbounds = NULL;
	}
	if (upb)
	{
		delete [] upb;
		upb = NULL;
	}
	if (lowb)
	{
		delete [] lowb;
		lowb = NULL;
	}
}

//template <class T> bool Optim::mod2DStartVal (int nparams, const double* start, const double* step, Matrix<T>* vals)
bool Optim::mod2DStartVal (int nparams, const double* start, const double* step, Matrix<double>* vals)
{
	if (nparams > 2)
		fprintf(stderr, "Too many parameters for Optim::mod2DStartVal\n");
	unsigned int npoints = 1;
	for (int i = 0; i < nparams; ++i)
	{
		if (step[i] > 0.0)
		{
			npoints *= Optim::numPoints(start[i], step[i]);
			if (npoints > vals->rown())
			{
				vals->~Matrix();
				vals->allocate(npoints, nparams);
			}
		}
		else if (step[i] == 0.0)
			continue;
		else
		{
			fprintf(stderr, "Invalid step size for parameter start position\n");
			return false;
		}
	}
	Optim::calc2DStart(nparams, npoints, start, step, vals);
	vals->setSubRow(npoints);
	vals->setSubCol(nparams);
	return true;
}

/*
template <class T> int Optim::guess2DStartVal (Matrix<T>* m, const T* guess, const T* p1lim, const T* p2lim, unsigned int npar1, unsigned int npar2)
{
	if (m->rown() != npar1*npar2 || m->coln() != 2)
	{
		fprintf(stderr, "Matrix must be of size npar1*npar2 X 2 in call to Optim::freq2DStartVal\n");
		return -1;
	}
	T p1inc = (p1lim[1] - p1lim[0]) / (npar1-1);
	T p2inc = (p2lim[1] - p2lim[0]) / (npar2-1);
	T p1start = guess[0];
	T p2start = guess[1];
	while (p1start - p1inc >= p1lim[0])
		p1start -= p1inc;
	while (p2start - p2inc >= p2lim[0])
		p2start -= p2inc;
	int j = 0;
	int k = 0;
	for (unsigned int i = 0; i < m->rown(); ++i)
	{
		if (!i % npar1)
			j = 0;
		if (!j % npar2)
			k = 0;
		(*m)[i][0] = p1start + j*p1inc;
		(*m)[i][1] = p2start + k*p2inc;
		++j;
		++k;
	}
	return 0;
}
*/

int Optim::guess2DStartVal (Matrix<double>* m, const double* guess , const double* p1lim, const double* p2lim, unsigned int npar1, unsigned int npar2)
{
	/*
	 * m: npar1*npar2 X 2 size matrix that holds values to be generated
	 * guess: a point in the matrix around which values will be generated
	 * p1lim: bounds for parameter 1 [lower, upper]
	 * p2lim: bounds for parameter 2 [lower, upper]
	 * npar1: number of unique values for parameter 1
	 * npar2: number of unique values for parameter 2
	 */
	if (m->rown() != npar1*npar2 || m->coln() != 2)
	{
		fprintf(stderr, "Matrix must be of size npar1*npar2 X 2 in call to Optim::freq2DStartVal\n");
		return -1;
	}
	double p1inc = npar1 == 1 ? 0 : (p1lim[1] - p1lim[0]) / (npar1-1);
	double p2inc = npar2 == 1 ? 0 : (p2lim[1] - p2lim[0]) / (npar2-1);
	double p1startval = guess[0];
	double p2startval = guess[1];
	if (npar1 > 1)
		while (p1startval - p1inc >= p1lim[0])
			p1startval -= p1inc;
	if (npar2 > 1)
		while (p2startval - p2inc >= p2lim[0])
			p2startval -= p2inc;
	int j = 0;
	int k = 0;
	double p1val = 0;
	double p2val = 0;
	for (unsigned int i = 0; i < m->rown(); ++i)
	{
		if (i > 0 && i % npar1 == 0)
			++j;
		if (i % npar2 == 0)
			k = 0;
		p1val = p1startval + j*p1inc;
		p2val = p2startval + k*p2inc;
		if (p1val > p1lim[1])
		{
			(*m)[i][0] = (p1lim[1] - (p1startval + (j-1)*p1inc) > guess[0] - p1lim[0]) ? p1val : p1lim[0];
		}
		else
			(*m)[i][0] = p1val;
		if (p2val > p2lim[1])
			(*m)[i][1] = (p2lim[1] - (p2startval + (k-1)*p2inc) > guess[1] - p2lim[0]) ? p2val : p2lim[0];
		else
			(*m)[i][1] = p2val;
		++k;
	}
	return 0;
}

Matrix<double> Optim::set2DStartVal (int nparams, const double* start, const double* step)
{
	Matrix<double> vals;
	if (nparams > 2)
	{
		fprintf(stderr, "Too many parameters for Optim::set2DStartVal\n");
		return vals;
	}
	unsigned int npoints = 1;
	for (int i = 0; i < nparams; ++i)
	{
		if (step[i] > 0.0)
			npoints *= Optim::numPoints(start[i], step[i]);
		else if (step[i] == 0.0)
			continue;
		else
		{
			fprintf(stderr, "Invalid step size for parameter start position\n");
			return vals;
		}
	}
	vals.allocate(npoints, nparams);
	Optim::calc2DStart(nparams, npoints, start, step, &vals);
	return vals;
}

template <class T> void Optim::calc2DStart(int nparams, unsigned int npoints, const double* start, const double* step, Matrix<T>* vals)
{
	unsigned int i;
	int j, k;
	int xbin = npoints/Optim::numPoints(start[0], step[0]);
	for (i = 0; i < npoints; ++i)
		for (j = 0; j < nparams; ++j)
			(*vals)[i][j] = start[j];

	for (i = 1; i < npoints; ++i)
	{
		if ( i % xbin == 0)
		{
			(*vals)[i][0] = (*vals)[i-1][0] + step[0];
			for (k = 1; k < nparams; ++k)
				(*vals)[i][j] = start[k];
		}
		else
		{
			(*vals)[i][0] = (*vals)[i-1][0];
			for (j = 1; j < nparams; ++j)
				(*vals)[i][j] = (*vals)[i-1][j] + step[j];
		}
	}
}

// finds number of points that lie in the interval [0,1] given a step size
int Optim::numPoints (double start, double step)
{

	int npts = 0;
	if (step == 0.0 || step > 1.0)
		return 1.0;
	else
	{
		double begin = 1 - start;
		if ( fmod(begin, step) == 0 )
			npts = begin/step + 1;
		else
			npts = ceil(begin/step);
	}
	return npts;
}

// runs the optimization at all starting locations
double Optim::runMultiOptim (double (*fn)(const double x[], const void*), void (*dfn)(const double x[], double y[], const void*),
	FunData* data, int optmethod, Matrix<double>* s)
{
	if (data->np > 2)
	{
		fprintf(stderr, "Too many parameters for Optim::runMultiOptim\n");
		negllh = -1.0;
		return negllh;
	}
	setExtraDat(data);
	double max = 1.0/0.0;
	Matrix<double>* startval = 0;
	if (s)
		startval = s;
	else
		startval = &start;
	size_t j = 0;
	size_t i = 0;
	static size_t param = 0;
	static size_t point = 0;
	point = startval->getSubRow() ? startval->getSubRow() : startval->rown();
	param = startval->getSubCol() ? startval->getSubCol() : startval->coln();
	double inpar [param];
	unsigned int enter_fail = _failcount;
	unsigned int exit_fail = optmethod ? (_failcount + point) : (_failcount + 2*point);
	while (i < point)
	{
		for (j = 0; j < param; ++j)
			inpar[j] = (*startval)[i][j];
		setParam(data->optp, inpar);
		if (!optmethod)
		{
			runBFGSB( fn, dfn ) ; // bfgs with analytic gradient
			if (_failcount >= enter_fail + point)
			{
				if (verb > 0)
					fprintf(stderr, "switching to numeric optimization for %i-parameter model\n", static_cast<int>(param));
				i = 0;
				optmethod = 1;
				continue;
			}
		}
		else
			runBFGSB( fn, NULL); // runs bfgs optimization with numerical gradient
		if (negllh < max)
		{
			max = negllh;
			for (j = 0; j < param; ++j)
				mlpar[j] = par[j];
		}
		++i;
	}
	if (_failcount < exit_fail)
		_failcount = 0;
	negllh = max;
	return negllh;
}


// returns _fail member
bool Optim::fail() const
{
	return _fail;
}
 // returns _nconditions
int Optim::conditions() const
{
	return _nconditions;
}

int Optim::getDim () const
{
	return dim;
}

double Optim::getNegllh () const
{
	return negllh;
}

double Optim::getParam (int index)
{
	if (par)
	{
		if (index < dim && index >= 0)
			return par[index];
		else
			fprintf(stderr, "index %i is out of parameter array range in Optim::getParam\n", index);
	}
	else
		fprintf(stderr, "parameter array not set in call to Optim::getParam\n");
	_fail = 1;
	return -1.0/0.0;
}

double Optim::getMLParam (int index)
{
	if (mlpar)
	{
		if (index < dim && index >= 0)
			return mlpar[index];
		else
			fprintf(stderr, "index %i is out of ML parameter array range in Optim::getMLParam\n", index);
	}
	else
		fprintf(stderr, "parameter array not set in call to Optim::getMLParam\n");
	_fail = 1;
	return -1.0/0.0;
}
