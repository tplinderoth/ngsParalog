// optim.h
// last edit: 3/3/2015

#ifndef OPTIM_H_
#define OPTIM_H_

#include "Matrix.h"

// STRUCTURES

// FunData stores information used by functions to optimize
struct FunData
{
	int np; // number parameters in model
	int optp; // number parameters to optimize
	void* sitedat; // sequencing data for a site
};

// CLASS DEFINITIONS
class Optim
{
	friend class Stats;
public:
	Optim (int noise = -1);
	~Optim ();
	bool setParN (int); // sets number of parameters to optimize
	bool setBoundCntrl (const double upperb [], int usize, const double lowerb [], int lsize, const int constraint [], int bsize); // sets boundary conditions
	bool setUpBound (const double u [], const int usize); // sets upper bound on params
	bool setLowBound (const double lo [], const int lsize); // sets lower bound on params
	bool setBound (const double upb [], int usize, const double lob [], int lsize); // sets upper and lower bounds
	bool setBoundCon (const int btype [], const int bsize); // sets the bounds constraint
	void setExtraDat (void *); // sets data used by likelihood function
	void setVerbose (int); // sets the verbosity of optimization
	bool setParam (const int npars, double* vals); // sets parameter starting values
	void setRandParam (const int npars); // sets random parameter starting values in [0,1]
	double runBFGSB (double (*fn)(const double x[], const void*), void (*dfn)(const double x[], double y[], const void*));
	void clearBounds ();
	static Matrix<double> set2DStartVal (int nparams, const double* start, const double* step); // creates matrix of optimization start points
	static int guess2DStartVal (Matrix<double>* m, const double* guess, const double* p1lim, const double* p2lim, unsigned int npar1, unsigned int npar2);
	//template <class T> static int guess2DStartVal (Matrix<T>* m, const T* guess, const T* xlim, const T* ylim, unsigned int npar1, unsigned int npar2);
	static bool mod2DStartVal (int nparams, const double* start, const double* step, Matrix<double>* vals); // modifies matrix of optimization start points
	double runMultiOptim (double (*fn)(const double x[], const void*), void (*dfn)(const double x[], double y[], const void*),
			FunData* data, int optmethod, Matrix<double>* s = 0); // runs bfgs at multiple starting points
	bool fail() const; // returns member _fail
	int conditions() const; // returns member _ncoditions
	int getDim() const; // returns member dim
	double getNegllh() const; // returns member negllh
	time_t getSeed();
	double getParam (int index); // returns parameter values
	double getMLParam (int index); // returns ML estimate of parameters
	Matrix<double> start; // 2D model start positions
private:
	static int numPoints (double start, double step); // finds number of points that lie in the interval [0,1] given a step size
	template<class T> void cpyConditions (int olddim, int newdim, T* oldarr); // copies existing optimization info when dimension is changed
	template <class T> static void calc2DStart(int nparams, unsigned int npoints, const double* start, const double* step, Matrix<T>* vals);
	void setMLPar();
	void seed();
	double* par; // parameter values
	double* mlpar; // maximum likelihood parameters
	double* upb; // upper bounds on the parameters
	double* lowb; // lower bounds on parameters
	int dim; // number of parameters to optimize
	int* nbounds; // boundary constraint option
	int verb; // amount of output by optimization
	void* xtradat; // data used in likelihood function
	int psize; // number of elements in par
	int lbsize; // number of elements in lowb
	int ubsize; // number of elements in upb
	int nbsize; // number of elements in nbounds
	bool _fail; // whether an error occurred in optimization
	int _nconditions; // number of conditions currently set (nparams, upbound, lowbound, params)
	double negllh; // negative log likelihood
	unsigned int _failcount; // number of optimization failures
	time_t _randset; // determines if srand has been called by optim object
};
#endif /* OPTIM_H_ */
