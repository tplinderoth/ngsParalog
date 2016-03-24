// Like.h

/*
 * for the parameter vector, p, of the alternative model p[0] = f, p[1] = m
 * for the null model p[0] = f
 */

#ifndef LIKE_H_
#define LIKE_H_

// CLASS DEFINITION

class Like
{
public:
	static double negLogfn (const double [], const void *);
	static void calcGradient (const double [], double [], const void *);
private:
	static double readProb (const double m, const double err, const int g1, const int g2, const int read);
	static void diffRead (const double err, const int g1, const int g2, const int read, double* dm, double* df);
	static double genoPrior (const double f, const int g1, const int g2);
	static double diffGenoPrior (const double f, const int g1, const int g2);
	static void kahanSum(double summand, double* total, double* comp);
	static void invalidRead ();
	static void invalidG1 ();
	static void invalidG2 ();
};
#endif /* LIKE_H_ */
