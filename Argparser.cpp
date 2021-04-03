/*
 * Argparser.cpp
 *
 *  Created on: Apr 19, 2016
 *      Author: tyler
 */

#include "Argparser.h"
#include "generalUtils.h"
#include <iomanip>
#include <cstring>

Argparser::Argparser ()
	: _task(0),
	  _minQ(20),
	  _Qoffset(33),
    _ploidy(2),
	  _minind(1),
	  _mincov(1),
	  _printml(0),
	  _numericGrad(1),
	  _verbose(-1),
	  _is(std::cin),
	  _os(std::cout.rdbuf()),
	  _fail(0)
{}

int Argparser::parseInput (const int c, char** v, const char* version)
{
// reads-in and checks input parameters

	int argPos = 1;
	int increment = 0;
	int rv = 0;

	if (c < 3)
	{
		if (c < 2)
		{
			maininfo(version);
			return 1;
		}
		else
			return (help(v[argPos],version));
	}

	_task = v[argPos];
	if ((rv = commandCheck(_task, v[2])) != 0)
		return rv;
	++argPos;

	while (argPos < c)
	{
		if (strcmp(v[argPos], "-infile") == 0)
		{
			_infile = (fexists(v[argPos + 1]) || strcmp(v[argPos + 1], "-") == 0) ? v[argPos + 1] : "";
			if (_infile.empty())
			{
				fprintf(stderr, "Couln't open input file %s\n",v[argPos+1]);
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-outfile") == 0)
			_outfile = v[argPos + 1];
		else if (strcmp(v[argPos], "-offsetQ") == 0)
		{
			_Qoffset = atof(v[argPos + 1]);
			if (_Qoffset < 0.0)
			{
				fprintf(stderr, "-offsetQ must be >= 0.0\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos],"-minQ") == 0)
		{
			_minQ = atof(v[argPos + 1]);
			if (_minQ < 0.0)
			{
				fprintf(stderr,"-minQ must be >= 0.0\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-minind") == 0)
		{
			_minind = atoi(v[argPos + 1]);
			if (_minind < 0)
			{
				fprintf(stderr, "-minind must be a positive integer\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-mincov") == 0)
		{
			_mincov = atoi(v[argPos + 1]);
			if (_mincov < 0)
			{
				fprintf(stderr, "-mincov must be a positive integer\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-printML") == 0)
		{
			_printml = atoi(v[argPos + 1]);
			if (_printml == 1 || _printml == 0)
			{}
			else
			{
				fprintf(stderr,"-printML must be enabled with 1 or disabled with 0\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-numericGrad") == 0)
		{
			_numericGrad = atoi(v[argPos + 1]);
			if (_numericGrad == 0 || _numericGrad == 1)
			{}
			else
			{
				fprintf(stderr,"-numericGrad 1 uses numerical gradient, -numericGrad 0 uses analytic gradient");
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-verbose") == 0)
    {
			_verbose = atoi(v[argPos+1]);
    }
    else if (strcmp(v[argPos], "-ploidy") == 0)
    {
      _ploidy = atoi(v[argPos+1]);
      if (_ploidy == 1 || _ploidy == 2)
      {}
      else
      {
        fprintf(stderr,"supported -ploidy are 1 or 2");
        return -1;
      }
    }
		else
		{
			fprintf(stderr, "Unknown argument: %s\n", v[argPos]);
			return -1;
		}
		argPos += 2 + increment;
		increment = 0;
	}

	return 0;
}

int Argparser::help (const char* arg, const char* version)
{
	if (strcmp(arg,"help") == 0 || strcmp(arg, "-help") == 0)
		maininfo(version);
	else if (strcmp(arg,"calcLR") == 0)
		calcLRinfo();
	else if (strcmp(arg,"findRegion") == 0)
		findRegionInfo();
	else
	{
		fprintf(stderr,"Unknown argument: %s\n",arg);
		return -1;
	}
	return 1;
}

int Argparser::commandCheck(const char* command, const char* arg1)
{
	if (strcmp(command,"calcLR") == 0)
	{
		if (strcmp(arg1,"-help") == 0)
		{
			calcLRinfo();
			return 1;
		}
		else
			return 0;
	}
	else if (strcmp(command,"findRegion") == 0)
		if (strcmp(arg1,"-help") == 0)
		{
			findRegionInfo();
			return 1;
		}
		else
			return 0;
	else
		fprintf(stderr,"Invalid command: %s\n",command);
	return -1;
}

int Argparser::setStreams (const std::string infile, const std::string outfile)
{
        // open input stream
        if (!infile.empty())
        {
                if (infile != "-")
                {
                        if (getFILE(_fin, infile.c_str(), "in")) // open the pileup file
                        {
                                _is.rdbuf(_fin.rdbuf()); // change input buffer from STDIN to fstream's buffer
                                std::cerr << "Reading input from " << infile << "\n";
                        }
                        else
                        {
                                std::cerr << "Problem opening input file " << infile << "\n";
                                _fail = 1;
                                return (_fail);
                        }
                }
                else
                std::cerr << "Reading from standard input\n";
        }
        else
        {
                std::cerr << "Problem parsing input stream\n";
                _fail=1;
                return(_fail);
        }

        // open output steam
        if (!outfile.empty())
        {
                if (getFILE(_fout, outfile.c_str(), "out"))
                {
                        _os.rdbuf(_fout.rdbuf()); // switch output stream's buffer to output file's buffer
                        std::cerr << "Dumping results to " << outfile << "\n";
                }
                else
                {
                        std::cerr << "Problem opening output file " << outfile << "\n";
                        _fail = 1;
                        return(_fail);
                }
        }
        else
                std::cerr << "Dumping results to standard output\n";

        return 0;
}

void Argparser::closeStreams ()
{
	if (_fin.is_open())
		_fin.close();
	if (_fout.is_open())
	{
		_fout.flush();
		_fout.close();
	}
}

unsigned int Argparser::ploidy () const
{
  return _ploidy;
}

double Argparser::minQ () const
{
	return _minQ;
}

double Argparser::offsetQ () const
{
	return _Qoffset;
}

unsigned int Argparser::minind () const
{
	return _minind;
}

unsigned int Argparser::mincov () const
{
	return _mincov;
}

std::string Argparser::infile_name () const
{
	return _infile;
}

std::string Argparser::outfile_name () const
{
	return _outfile;
}

std::istream& Argparser::input ()
{
	return _is;
}

std::ostream& Argparser::output ()
{
	return _os;
}

char* Argparser::task () const
{
	return _task;
}

int Argparser::printML() const
{
	return _printml;
}

int Argparser::isNumericGrad () const
{
	return _numericGrad;
}

int Argparser::verblevel () const
{
	return _verbose;
}

int Argparser::fail () const
{
	return _fail;
}

void Argparser::maininfo (const char* v)
{
	int w = 12;
	std::cerr << "\nngsParalog\nversion " << v << "\n\nUsage: ngsParalog [command] [arguments]\n"
	<< "\nCommands:\n"
	<< "\n" << std::setw(w) << std::left << "calcLR" << "Calculate per site likelihood ratio of duplication"
	<< "\n" << std::setw(w) << std::left << "findRegion" << "Find the coordinates of duplicated genomic regions"
	<< "\n\n";
}

void Argparser::calcLRinfo ()
{
	int w = 12;
	std::cerr << "\nUsage: ngsParalog calcLR [arguments]\n"
	<< "\nArguments:\n"
	<< "\n" << std::setw(w) << std::left << "-infile" << std::setw(w) << "<string>" << "Pileup format file of reads and quality scores; specifying '-' will read from STDIN [" << _infile << "]"
	<< "\n" << std::setw(w) << std::left << "-outfile" << std::setw(w) << "<string>" << "Name of output file, if not provided results printed to STDOUT [" << _outfile << "]"
	<< "\n" << std::setw(w) << std::left << "-offsetQ" << std::setw(w) << "<float>" << "Quality score offset [" << _Qoffset << "]"
	<< "\n" << std::setw(w) << std::left << "-minQ" << std::setw(w) << "<float>" << "Minimum base quality score [" << _minQ << "]"
	<< "\n" << std::setw(w) << std::left << "-minind" << std::setw(w) << "<int>" << "Minimum number of covered individuals to retain site [" << _minind << "]"
	<< "\n" << std::setw(w) << std::left << "-mincov" << std::setw(w) << "<int>" << "Minimum number of reads for an individual to be considered 'covered' [" << _mincov << "]"
	<< "\n\nOutput by field:"
	<< "\n(1) sequence ID"
	<< "\n(2) position in sequence (1-base indexed)"
	<< "\n(3) -log null likelihood"
	<< "\n(4) -log alternative likelihood"
	<< "\n(5) likelihood ratio of mismapped reads"
	<< "\n\n";
}

void Argparser::findRegionInfo ()
{
	//int w = 12;
	std::cerr << "\nUsage: ngsParalog findRegion [arguments]\n"
	<< "\nArguments:\n"
	<< "\n\nOutput:\n"
	<< "\n";
}
