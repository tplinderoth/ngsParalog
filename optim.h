// optim.h
// last edit: 3/3/2015

#ifndef OPTIM_H_
#define OPTIM_H_

#include "Matrix.h"

// CLASS DEFINITIONS
class Optim
{
public:
	Optim (int noise = -1);
	~Optim ();

	// set member functions
	bool setParN (int n, bool init = true); // sets number of parameters to optimize
	bool setBoundCntrl (const double upperb [], int usize, const double lowerb [], int lsize, const int constraint [], int bsize); // sets boundary conditions
	bool setUpBound (const double u [], const int usize); // sets upper bound on params
	bool setLowBound (const double lo [], const int lsize); // sets lower bound on params
	bool setBound (const double upb [], int usize, const double lob [], int lsize); // sets upper and lower bounds
	bool setBoundCon (const int btype [], const int bsize); // sets the bounds constraint
	void setData (void *); // sets data used by likelihood function
	void setVerbose (int); // sets the verbosity of optimization
	bool setParam (const int npars, double* vals = NULL); // sets parameter starting values
	void setRandParam (const int npars); // sets random parameter starting values in [0,1]
	double& setllh (); // assign value to _llh member
	void setmlparam (int index, double val); // assign value to _mlpar
	int initStartMatrix (int nparams, const double* s, const double* step); // set points in start member

	// return member functions
	double* lowbounds () const; // return pointer to member _lowb
	double* upbounds () const; // return pointer to member _upb
	int* numbounds () const; // return point to member _nbounds
	int verblevel () const; // return member _verb
	int* fail(); // point to member _fail
	int isfail() const; // returns _fail member
	double llh() const; // returns member _llh
	int conditions() const; // returns member _ncoditions
	time_t getSeed();
	double getParam (int index); // returns parameter values
	double mlparam (int index); // returns ML estimate of parameters
	int getDim() const; // returns member _dim

	// operational functions
	void clearBounds ();
	static Matrix<double> genStartMatrix (int nparams, const double* start, const double* step); // creates matrix of optimization start points
	static int guess2DStartVal (Matrix<double>* m, const double* guess, const double* p1lim, const double* p2lim, unsigned int npar1, unsigned int npar2);
	//template <class T> static int guess2DStartVal (Matrix<T>* m, const T* guess, const T* xlim, const T* ylim, unsigned int npar1, unsigned int npar2);
	static bool mod2DStartVal (int nparams, const double* start, const double* step, Matrix<double>* vals); // modifies matrix of optimization start points
	double multiOptim (double (*fn)(const double x[], const void*), void (*dfn)(const double x[], double y[], const void*),  Matrix<double>* s = NULL);

	// public member variables
	Matrix<double> start; // 2D model start positions
	double* par; // parameter values
	void* data; // data used in likelihood function

private:
	// functions
	static int numPoints (double start, double step); // finds number of points that lie in the interval [0,1] given a step size
	template<class T> void cpyConditions (int olddim, int newdim, T* oldarr); // copies existing optimization info when dimension is changed
	template <class T> static void calc2DStart(int nparams, unsigned int npoints, const double* start, const double* step, Matrix<T>* vals);
	void initmlparam(int n, double* vals = NULL);
	void seed();

	//private member variables
	double* _mlpar; // maximum likelihood parameters
	double* _upb; // upper bounds on the parameters
	double* _lowb; // lower bounds on parameters
	int _dim; // number of parameters to optimize
	int* _nbounds; // boundary constraint option
	int _verb; // amount of output by optimization
	int _fail; // whether an error occurred in optimization
	int _nconditions; // number of conditions currently set (nparams, upbound, lowbound, params)
	double _llh; // negative log likelihood
	time_t _randset; // determines if srand has been called by optim object
};
#endif /* OPTIM_H_ */
