/*
 * Argparser.h
 *
 *  Created on: Apr 19, 2016
 *      Author: tyler
 */

#ifndef ARGPARSER_H_
#define ARGPARSER_H_

#include <fstream>
#include <string>
#include <cstdlib>

class Argparser
{
public:
        /* FUNCTIONS */
        Argparser ();
        int parseInput (const int c, char** v, const char* version); /* parse command line arguments to set data members*/
        int setStreams (std::string infile, std::string outfile); /* sets IO streams */
        void closeStreams (); /* close IO streams */
        unsigned int ploidy () const;
        double minQ () const; /* return _minQ */
        double offsetQ () const; /* return _Qoffset */
        unsigned int minind () const; /* return minind */
        unsigned int mincov () const; /* return _mincov */
        std::string infile_name () const; /* returns _infile */
        std::string outfile_name () const; /* returns _outfile */
        std::istream& input (); /* returns &_is */
        std::ostream& output (); /* returns _os */
        char* task () const; /* return _task */
        int printML () const; /* return _printml */
        int isNumericGrad () const; /* return _numericGrad */
        int verblevel () const; /* return _verbose */
        int fail () const; /* return _fail */
private:
        /* FUNCTIONS */
        int help (const char* arg, const char* version); /* calls subroutine help */
        void maininfo (const char* v); /* print main help */
        void calcLRinfo (); /* print help information for 'calcLR' subroutine */
        void findRegionInfo (); /* print help information for 'findRegion' subroutine */
        int commandCheck(const char* command, const char* arg1); /* checks for valid commands */
        /* MEMBER VARIABLES */
        char* _task; /* 0 - calcLR, 1 - findRegion */
        std::string _infile; /* name of input pileup file */
        std::string _outfile; /* name of output file */
        double _minQ; /* min base quality score */
        double _Qoffset; /* min possible ASCII decimal value used to encode quality scores */
        unsigned int _ploidy; /* ploidy of organism */
        int _minind; /* min number of covered individuals */
        int _mincov; /* min number of reads for an individual to be considered "covered" */
        int _printml; /* controls what parameter ML estimates are outputted */
        int _numericGrad; /* 1 - use numeric gradient, 0 - use analytic gradient */
        int _verbose; /* controls amount optimization output */
        /* IO streams */
        std::fstream _fin;
        std::fstream _fout;
        std::istream& _is;
        std::ostream _os;
        /* status flags */
        int _fail; /* flag to denote program errors */
};

#endif /* ARGPARSER_H_ */
