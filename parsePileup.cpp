// parsePileup.cpp

#include <cstdio>
#include <cstdlib>
#include "generalUtils.h"
#include "parsePileup.h"

// SiteData structure constructor
SiteData::SiteData ()
	: depth(0)
{}

// Pileup class constructor
Pileup::Pileup ()
	: name("none"),
	  pos(0),
	  encode(33),
	  minQ(13),
	  _fail(0),
	  _nind(0),
	  _depthReserve(20),
	  _numalt(0),
	  _numref(0),
	  _alleles(),
	  _refallele('\0')
{}

// Pileup::setQualCode assigns quality encoding
void Pileup::setQualCode (int code)
{
	if (code < 0)
	{
		fprintf(stderr, "Invalid quality score encoding in Pileup::getCode\n");
		_fail = 1;
	}
	encode = code;
}

// Pileup::setMinQ assigns minimum quality score for reads to be considered
void Pileup::setMinQ (int q)
{
	if (q < 0)
	{
		fprintf(stderr, "Cannot set minimum quality score to < 0 in Pileup::setMinQ\n");
		_fail = 1;
	}
	minQ = q;
}

int Pileup::getMinQ ()
{
 return minQ;
}

int Pileup::getQualCode ()
{
	return encode;
}

// Pileup::getSeqDat extracts information from pileup line
int Pileup::getSeqDat (const std::string& pile)
{
	if (pile.empty())
	{
		fprintf(stderr, "No pileup line supplied to Pileup::getSeqDat\n");
		_fail = 1;
	}

	std::vector<std::string> pvec = split(pile, '\t');

	if (pvec.size() == 1)
	{
		fprintf(stderr, "Piluep line couldn't be split by Pileup::getSeqDat; check delimiter\n");
		_fail = 1;
	}

	// assign name
	if (!pvec[0].empty())
		name = pvec[0];
	else
	{
		fprintf(stderr, "Attempt to assign empty container to name by Pileup::getSeqDat\n");
		_fail = 1;
	}

	// assign position
	if (!pvec[1].empty())
		pos = atoi(pvec[1].c_str());
	else
	{
		fprintf(stderr, "Attempt to assign empty container to position by Pileup::getSeqDat\n");
		_fail = 1;
	}

	// get reference allele identity
	if (!pvec[2].empty())
		_refallele = toupper(pvec[2][0]);
	else
	{
		fprintf(stderr, "Attempt to assign empty container as reference allele by Pileup::getSeqDat\n");
		_fail = 1;
	}

	// get number of individuals
	if (seqdat.empty())
	{
		_nind = ExtractIndN(pile);
		if (_nind < 1)
		{
			std::cerr << "Number of individuals < 1\n";
			return -1;
		}
		initializeSeqdat(_nind);
	}

	// reset allele counts
	if (_numref + _numalt > 0)
	{
		_numref = 0;
		_numalt = 0;
		for (int i = 0; i < 5; ++i)
			_alleles[i] = 0;
	}

	// assign reads and quality scores
	size_t ind = 0;
	std::string reads;
	std::string qscores;
	std::vector<std::string>::const_iterator iter = pvec.begin() + 3;

	if ((*iter).empty())
	{
		fprintf(stderr, "Sequencing data required by Pileup::getSeqDat missing\n");
		_fail = 1;
	}

	for (iter = pvec.begin() + 3; iter != pvec.end(); ++iter)
	{
		++ind;
		seqdat[ind-1].depth = 0;
		if ( (*iter)[0] == '0' ) // no data for individual
		{
			while ((iter+1) != pvec.end() && !isdigit((*(iter+1))[0]))
				++iter;
			if ((iter+1) == pvec.end())
				break;
		}
		else
		{
			reads = *(++iter);
			while ( *((iter)->end()-1) == '^' ) // ascii "space" is mapping quality
			{
				reads += ' ';
				reads += *(++iter);
			}
			qscores = *(++iter);
			getReadDat(&reads, &qscores, ind-1);
		}
	}
	return 0;
}

// assigns read and quality scores for an individual
void Pileup::getReadDat (std::string* reads, std::string* qual, size_t ind)
{
	char r;
	int phredq = 0;
	unsigned int index = 0;
	unsigned int indsize = 0;
	size_t read_num = 0;
	std::string::iterator q = qual->begin();

	for(std::string::const_iterator r_iter = reads->begin(); r_iter != reads->end(); ++r_iter) // go over all reads
	{
		r = toupper(*r_iter);
		if (r == '.' || r == ',') // ref
		{
			phredq = *q - encode;
			if (phredq >= minQ)
			{
				newRead(ind, read_num, 0, phredq);
				++read_num;
				++_numref;
				switch(_refallele)
				{
					case 'A' :
						++_alleles[0];
						break;
					case 'C' :
						++_alleles[1];
						break;
					case 'G' :
						++_alleles[2];
						break;
					case 'T' :
						++_alleles[3];
						break;
					default :
						fprintf(stderr, "warning: unrecognized reference allele '%c' at %s position %d\n", _refallele, name.c_str(), pos);
				}
			}
			++q;
			++index;
		}
		else if (r == 'A' || r == 'C' || r == 'G' || r == 'T') // alt
		{
			phredq = *q - encode;
			if (phredq >= minQ)
			{
				newRead(ind, read_num, 1, phredq);
				++read_num;
				++_numalt;
				switch(r)
				{
					case 'A' :
						++_alleles[0];
						break;
					case 'C' :
						++_alleles[1];
						break;
					case 'G' :
						++_alleles[2];
						break;
					case 'T' :
						++_alleles[3];
						break;
				}
			}
			++q;
			++index;
		}
		else if (r == 'N' || r == '*' || r == '<' || r == '>') // missing read or gaps
		{
			++q;
			++index;
		}
		else if (r == '+' || r == '-') // indel
		{
			indsize = indelSize(reads, (index+1));
			r_iter += indsize;
			index += indsize + 1;
			++_alleles[4];
		}
		else if (r == '^') // beginning of read
		{
			++r_iter;
			index += 2;
		}
		else if (r == '$') // end of read
			++index;
		else
		{
			fprintf(stderr, "Unrecognized symbol '%c' in pileup line by Pileup::getReadDat\n", r);
			std::cerr << "bad site: " << name << " " << pos << "\n" << *reads << " " << *qual << "\n";
			_fail = 1;
			//exit(0); // exit point for debugging
		}
	}
	seqdat[ind].depth = read_num;
}

// Pileup::indelSize gets the size of an indel + the number of digits comprising the size
unsigned int Pileup::indelSize (std::string* s, unsigned int start)
{
        unsigned int len = 0;
        std::string::const_iterator it = s->begin() + start;
        while ( isdigit(*it) && it != s->end())
        {
                ++len;
                ++it;
        }

        std::string indsize = s->substr (start,len);
        return ( atoi(indsize.c_str()) + indsize.length() );
}

// Pileup::fail returns value of _fail member
int Pileup::fail ()
{
	return _fail;
}

size_t Pileup::ExtractIndN (const std::string& line)
{
	size_t nind = 0;
	int inc = 0;
	if (!line.empty())
	{
		std::vector<std::string> ptoke = split(line, '\t');
		for (std::vector<std::string>::const_iterator iter = ptoke.begin() + 3; iter != ptoke.end(); ++iter)
		{
			++nind;
			inc = 0;
			if ( (*iter)[0] == '0' ) // no data for individual
			{
				while ((iter+1) != ptoke.end() && !isdigit((*(iter+1))[0]))
					++iter;
				if ((iter+1) == ptoke.end())
					break;
			}
            		else
            		{
            			while ( *((iter+(1+inc))->end()-1) == '^' ) // ascii "space" is mapping quality
            			{
            				++inc;
            			}
            			iter += 2 + inc;
            		}
		}
	}
	else
	{
		fprintf(stderr, "No sequencing data in string provided to Pileup::ExtractIndN\n");
		_fail = 1;
	}
	return nind;
}

void Pileup::setIndN (size_t n)
{
	if (n > 0)
		_nind = n;
	else
		fprintf(stderr, "Attempt to set nonpositive number of individuals in Pileup::setIndN\n");
}

void Pileup::initializeSeqdat (size_t n)
{
	if (n >= 0)
	{
		size_t start_size = _depthReserve;
		seqdat.resize(n);
		for (size_t i = 0; i < n; ++i)
		{
			seqdat[i].rdat.resize(start_size);
			seqdat[i].depth = start_size;
		}
	}
}

void Pileup::newRead (const size_t ind, const size_t read_num, int read, int quality)
{
	static const float factor = 0.5;
	if (quality < 0)
	{
		fprintf(stderr, "ERROR: Negative quality score '%i' --> check quality score encoding\n", quality);
		_fail = 1;
	}
	if (read_num >= seqdat[ind].rdat.size())
	{
		if (read_num < seqdat[ind].rdat.capacity())
		{
			seqdat[ind].rdat.push_back(std::make_pair(read, quality));
			return;
		}
		else
			seqdat[ind].rdat.resize( seqdat[ind].rdat.size() + factor * seqdat[ind].rdat.size() );
	}
	seqdat[ind].rdat[read_num].first = read;
	seqdat[ind].rdat[read_num].second = quality;
}

std::string Pileup::seqName ()
{
	return name;
}

unsigned int Pileup::position ()
{
	return pos;
}

void Pileup::setDepthReserve (size_t depth)
{
	_depthReserve = depth;
}

size_t Pileup::nInd ()
{
	return _nind;
}

size_t Pileup::altcount ()
{
	return _numalt;
}

size_t Pileup::refcount ()
{
	return _numref;
}

double Pileup::altfreq ()
{
	return  static_cast <double> (_numalt) / (_numref + _numalt);
}

size_t Pileup::siteDepth ()
{
	return _numref + _numalt;
}

unsigned int Pileup::alleleCount (char allele)
{
	unsigned int n = 0;
	switch (allele)
	{
		case 'A' :
		case 'a' :
			n = _alleles[0];
			break;
		case 'C' :
		case 'c' :
			n = _alleles[1];
			break;
		case 'G' :
		case 'g' :
			n = _alleles[2];
			break;
		case 'T' :
		case 't' :
			n = _alleles[3];
			break;
		case 'I' :
		case 'i' :
			n = _alleles[4];
			break;
		default :
			fprintf(stderr, "Unrecognized base '%c' in call to Pileup::alleleCount\n", allele);
	}
	return n;
}

char Pileup::refAllele ()
{
	return _refallele;
}
