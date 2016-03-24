/*
 * stats.h
 *
 */

#ifndef STATS_H_
#define STATS_H_

#include "optim.h"

// CLASS DEFINITIONS

class Stats
{
public:
	Stats();
	~Stats();
	double optimLR ( double (*fun)(const double x[], const void*), void (*dfun)(const double x[], double y[], const void*),
			FunData* likedat, Optim* altopt, Optim* nullopt, int gradMeth, bool logscale, const int* reduced_ind, double increment,
			double decrease, int attempt, int tries, double lr_threshold, Matrix<double>* alt_start = 0, Matrix<double>* null_start = 0, int verb = -1);
/*
	* optimLR:
	* fun <- likelihood function
	* dfun <- function that calculates the gradient of the likelihood function
	* lrfun <- function that calculates the likelihood ratio
	* gradMeth <- variable to indicate whether to use numeric or analytic gradient
	* logscale <- variable to indicate whether likelihood returned by fun is in -log likelihood scale
	* reduced_ind <- an array that stores the INDEX of null model parameters in the full model parameter array
	* increment <- how much to increment parameter value in starting position matrix
	* decrease <- factor to decrease the increment amount by in subsequent attempts to optimize (produces more start points)
	* attempt <- how many attempts to optimize with increment/decrease method before generating random start points
	* tries <- number of failed optimization attempts before giving up and moving on
	* lr_threshold <- threshold for rejecting LR, protects against optimization failure
	* alt_start <- optional alternative model start points
	* null_start <- optional null model start points
	* verb <- controls the amount of output
*/
	double multiLR (double (*fun)(const double x[], const void*), void (*dfun)(const double x[], double y[], const void*),
			FunData* likedat, Optim* altopt, Optim* nullopt, int gradMeth, bool logscale,
			Matrix<double>* alt_start = 0, Matrix<double>* null_start = 0);
	static double calcLR (const double null, const double alt, bool islog = true); // if islog = 1 likelihoods are in -log scale, otherwise islog = 0
	bool fail (); // returns good_calc
	size_t getFullStepSize ();
	size_t getReducedStepSize ();
	size_t getFullStartSize ();
	size_t getReducedStartSize ();
	void seed();
	time_t getSeed();
private:
	void setFullStepSize (int n);
	void setReducedStepSize (int n);
	void setFullStartSize (int n);
	void setReducedStartSize (int n);
	void optFailMsg (FunData* data, const int thresh);
	bool _fail; // flags indicates "true" if calculation was successful or "false" if error
	double* _step_full;
	double* _step_reduced;
	double* _start_full;
	double* _start_reduced;
	int _step_full_size;
	int _step_reduced_size;
	int _start_full_size;
	int _start_reduced_size;
	double _stat;
	time_t _seed;
};

#endif /* STATS_H_ */
