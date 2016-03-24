// generalUtils.cpp
// version 0.0.8; 15 December 2014

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <sstream>
#include "generalUtils.h"

// getCString converts string to C-Style pointer string
char * getCString (std::string s)
{
 	int slen = s.length();
	char * cstring = new char[slen + 1];
	s.copy(cstring, slen, 0);
	cstring[ slen ] = '\0';

	return cstring;
}

// getFILE is a wrapper for getting files
bool getFILE(std::fstream &fp, const char* fname, const char* mode)
{
	int writeFile = 0;
	if (strcmp(mode, "out") == 0)
	{
		writeFile = 1;
		if(writeFile && fexists(fname))
		{
			fprintf(stderr,"File already exists: %s\n",fname);
			return false;
		}

		fp.open(fname, std::ios::out);
	}
	else if (strcmp(mode, "app") == 0)
		fp.open(fname, std::ios::app);
	else if (strcmp(mode, "in") == 0)
		fp.open(fname, std::ios::in);

	if( !fp )
	{
		fprintf(stderr,"Error opening FILE handle for file: %s\n",fname);
		fp.close();
		return false;
	}
	return true;
}

// fexists finds out if a file exists
int fexists(const char* str)
{
	struct stat buffer;
 	return (stat(str, &buffer )==0 );
}

// readChunk chunks in a file line by line (not efficient since it doesn't actually read a chunk)
bool readChunk (std::vector<std::string>& datavec, unsigned int* chunk, int* end, std::istream& is)
{
	unsigned int i = 0;
	std::string(line);
	size_t cap = 0;

	if (!datavec.empty())
	{
		cap = datavec.capacity();
		datavec.clear();
		if (datavec.capacity() < cap)
			datavec.reserve(cap);
	}
	while (i < *chunk)
	{
		if (getline(is, line))
		{
			datavec.push_back(line);
			++i;
		}
		else if (is.eof())
		{
			datavec.resize(i); // remove empty elements
			*end = 1;
			*chunk = i; // number of lines processed upon returning
			return(true);
		}
		else
		{
			fprintf(stderr, "Failed reading stream in readChunk function\n");
			return(false);
		}
	}
	return(true);
}

// split splits a string based on a delimiter
std::vector<std::string> split (const std::string& s, char delim)
{
        std::vector<std::string> elems;
        std::stringstream ss(s);
        std::string sholder;
        while (std::getline(ss, sholder, delim))
        {
                if (!sholder.empty())
                        elems.push_back(sholder);
        }
        return elems;
}

// draws random, bounded, decimal number
double decimalUnifBound (double min, double max )
{
        return rand() / (static_cast<double>(RAND_MAX) + 1) * (max - min) + min;
}
