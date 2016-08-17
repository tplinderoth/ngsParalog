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
	  data(0),
	  _mlpar(0),
	  _upb(0),
	  _lowb(0),
	  _dim(0),
	  _nbounds(0),
	  _fail(false),
	  _nconditions(0),
	  _llh(0.0),
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
	if (_mlpar)
		delete [] _mlpar;

	//clear bounds
	if (_nbounds)
	{
		delete [] _nbounds;
		--_nconditions;
	}
	if (_upb)
	{
		delete [] _upb;
		--_nconditions;
	}
	if (_lowb)
	{
		delete [] _lowb;
		--_nconditions;
	}
}


template< class T > void Optim::cpyConditions (int olddim, int newdim, T* oldarr)
{
// copies existing optimization info when dimension is changed

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

bool Optim::setParN (int n, bool init)
{
// gets number of parameters to optimize

	if (n > 0)
	{
		if (_dim)
		{
			if (par)
				cpyConditions(_dim, n, par);
			if (_lowb)
				cpyConditions(_dim, n, _lowb);
			if (_upb)
				cpyConditions(_dim, n, _upb);
			if (_nbounds)
				cpyConditions(_dim, n, _nbounds);
			if (_mlpar)
				cpyConditions(_dim, n, _mlpar);
		}
		_dim = n;
		if (!par && init)
			setParam(_dim);
	}
	else
	{
		fprintf(stderr, "Invalid number of parameters in Optim::setParN\n");
		_dim = 0;
		return false;
	}
	return true;
}


bool Optim::setBoundCntrl (const double upperb [], int usize, const double lowerb [], int lsize, const int constraint [], int bsize)
{
// sets boundary conditions for optimization

	if (!setBoundCon(constraint, bsize))
		return false;
	if (!setBound(upperb, usize, lowerb, lsize))
			return false;
	return true;
}

bool Optim::setUpBound (const double u [], const int usize)
{
// sets the upper bound on parameters to optimize

	if (usize == _dim)
	{
		if (!_upb)
		{
			_upb = new double [_dim];
			++_nconditions;
		}
		for (int i = 0; i < _dim; ++i)
			_upb[i] = u[i];
	}
	else
	{
		fprintf(stderr, "Number of bounds doesn't match Optim::dim in Optim::setUpBound\n");
		_upb = NULL;
		return false;
	}
	return true;
}

bool Optim::setLowBound (const double lo [], const int lsize)
{
// sets the lower bound on parameters to optimize

	if (lsize == _dim)
	{
		if (!_lowb)
		{
			_lowb = new double [_dim];
			++_nconditions;
		}
		for (int i = 0; i < _dim; ++i)
			_lowb[i] = lo[i];

	}
	else
	{
		fprintf(stderr, "Number of bounds doesn't match Optim::dim in Optim::setLowBound\n");
		_lowb = NULL;
		return false;
	}
	return true;
}

void Optim::initmlparam(int n, double* vals)
{
	if (_dim > 0)
		if (!_mlpar)
        {
			_mlpar = new double[n];
            for (int i = 0; i < n; ++i)
            	_mlpar[i] = vals ? vals[i] : 0.0;
            ++_nconditions;
        }
        else
        	fprintf(stderr, "Optim::_mlpar already initiated in call to Optim::initmlparam()\n");
	else
		fprintf(stderr, "warning: Optim::dim not set in call to Optim::initmlparam\n");
}

bool Optim::setBound (const double upb [], int usize, const double lob [], int lsize)
{
// sets the upper and lower boundaries on parameters to optimize

	if (!setUpBound(upb, usize)) // set upper bound
		return false;
	if (!setLowBound(lob, lsize)) // set lower bound
		return false;
	return true;
}

bool Optim::setBoundCon (const int btype [], const int bsize)
{
// sets the number of boundaries for each parameter

	if (bsize == _dim)
	{
		if (!_nbounds)
		{
			_nbounds = new int [_dim];
			++_nconditions;
		}
		for (int i = 0; i < _dim; ++i)
			_nbounds[i] = btype[i];
	}
	else
	{
		fprintf(stderr, "Number of bound constraint types doesn't match number of parameters in Optim::setBoundCon\n");
		_nbounds = NULL;
		return false;
	}
	return true;
}


void Optim::setVerbose (int v)
{
// sets the amount of output during optimization
	_verb = v;
}

void Optim::setData (void* x)
{
// sets data used by likelihood function
	if (x)
		data = x;
	else
	{
		fprintf(stderr, "Assignment of NULL pointer in Optim::setExtraDat\n");
		data = NULL;
	}
}

bool Optim::setParam (const int npars, double* vals)
{
// sets user defined starting parameter values
	if (npars == _dim)
	{
		if (!par)
		{
			par = new double [_dim];
			++_nconditions;
		}
		for (int i = 0; i < _dim; ++i)
			par[i] = vals ? vals[i] : 0;
		if (!_mlpar)
			initmlparam(_dim);
	}
	else
	{
		fprintf(stderr, "Number of parameters supplied to Optim::setParam differs from Optim::dim\n");
		par = NULL;
		return false;
	}
	return true;
}

void Optim::setRandParam (const int npars)
{
// sets random starting parameter values in [0,1]

	if (!_randset)
		seed();
	if (npars == _dim)
	{
		if (!par)
		{
			par = new double [_dim];
			++_nconditions;
		}
		for (int i = 0; i < _dim; ++i)
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

void Optim::clearBounds ()
{
// frees memory allocated for bounding

	if (_nbounds)
	{
		delete [] _nbounds;
		_nbounds = NULL;
	}
	if (_upb)
	{
		delete [] _upb;
		_upb = NULL;
	}
	if (_lowb)
	{
		delete [] _lowb;
		_lowb = NULL;
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

Matrix<double> Optim::genStartMatrix (int nparams, const double* start, const double* step)
{
	Matrix<double> vals;
	if (nparams > 2)
	{
		fprintf(stderr, "Too many parameters for Optim::genStartMatrix\n");
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

int Optim::initStartMatrix (int nparams, const double* s, const double* step, int extrapoints)
{
	if (start.rown() != 0 || start.coln() != 0)
	{
		fprintf(stderr, "Optim start member already initialized\n");
		return 1;
	}
	if (nparams > 2)
	{
		fprintf(stderr, "Too many parameters for Optim::initStartMatrix\n");
		_fail = 1;
		return 1;
	}
	unsigned int npoints = 1;
	for (int i = 0; i < nparams; ++i)
	{
		if (step[i] > 0.0)
			npoints *= Optim::numPoints(s[i], step[i]);
		else if (step[i] == 0.0)
			continue;
		else
		{
			fprintf(stderr, "Invalid step size for parameter start position\n");
			_fail = 1;
			return 1;
		}
	}
	start.allocate(npoints+extrapoints, nparams); // add 1 row to place in null condition, add and isalt parameter +1 to npoints parameter
	Optim::calc2DStart(nparams, npoints, s, step, &start);
	++_nconditions;
	return 0;
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


int Optim::numPoints (double start, double step)
{
// finds number of points that lie in the interval [0,1] given a step size

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

double Optim::multiOptim (double (*fn)(const double x[], const void*), void (*dfn)(const double x[], double y[], const void*), Matrix<double>* s)
{
/*
 * run optimization at all start points using L-BFGS-B algorithm
 * fn is function to be minimized
 * dfn is the funtion to calculate the gradient of fn, dfn=NULL uses numeric gradient
 */
	_fail = 0;
	if (_dim > 2)
		fprintf(stderr, "Too many parameters for Optim::multiOptim\n");
	if (data == NULL)
		fprintf(stderr, "Optim data member not set in call to Optim::multiOptim\n");
	if (_fail)
	{
		_llh = 0.0;
		_fail = 1;
		return _llh;
	}
	double max = 1.0/0.0;
	int success = 0;
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
	while (i < point)
	{
		for (j = 0; j < param; j++)
			inpar[j] = (*startval)[i][j];
		setParam(_dim, inpar);
		_llh = findmax_bfgs(_dim, par, this, fn, dfn, _lowb, _upb, _nbounds, _verb, &_fail);
		if (_llh < max && !_fail)
		{
			++success;
			max = _llh;
			for (j = 0; j < param; ++j)
				_mlpar[j] = par[j];
		}
		_fail = 0;
		++i;
	}
	_fail = success ? 0 : 1;
	_llh = max;
	return _llh;
}

int Optim::isfail() const
{
	return _fail;
}

int* Optim::fail()
{
	return &_fail;
}

int Optim::conditions() const
{
	return _nconditions;
}

int Optim::getDim () const
{
	return _dim;
}

double Optim::llh () const
{
	return _llh;
}

double Optim::getParam (int index)
{
	if (par)
	{
		if (index < _dim && index >= 0)
			return par[index];
		else
			fprintf(stderr, "index %i is out of parameter array range in Optim::getParam\n", index);
	}
	else
		fprintf(stderr, "parameter array not set in call to Optim::getParam\n");
	_fail = 1;
	return -1.0/0.0;
}

double Optim::mlparam (int index)
{
	if (_mlpar)
	{
		if (index < _dim && index >= 0)
			return _mlpar[index];
		else
			fprintf(stderr, "index %i is out of ML parameter array range in Optim::mlparam\n", index);
	}
	else
		fprintf(stderr, "parameter array not set in call to Optim::mlparam\n");
	_fail = 1;
	return 0.0;
}

void Optim::setmlparam (int index, double val)
{
	if (_mlpar)
	{
		if (index < _dim && index >= 0)
			_mlpar[index] = val;
		else
		{
			fprintf(stderr, "index %i is out of ML parameter array range in Optim::mlparam\n", index);
			_fail = 1;
		}
	}
	else
	{
		fprintf(stderr, "_mlparam not initiated in call to Optim::mlparam\n");
		_fail = 1;
	}
}

double& Optim::setllh ()
{
	return _llh;
}

double* Optim::lowbounds () const
{
	return _lowb;
}

double* Optim::upbounds () const
{
	return _upb;
}

int Optim::verblevel () const
{
	return _verb;
}

int* Optim::numbounds () const
{
	return _nbounds;
}
