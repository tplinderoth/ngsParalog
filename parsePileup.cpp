// parsePileup.cpp

#include <cstdio>
#include <algorithm> // sort
#include <math.h> // pow
#include <map>
#include <cstring> // strcpy
#include "parsePileup.h"
#include "generalUtils.h"


// SiteData structure constructor
SiteData::SiteData ()
	: depth(0)
{}

unsigned int SiteData::cov (char allele) const
{
	if (allele)
	{
		switch(allele)
		{
			case 'A' :
			case 'a' :
				return allecount[0];
			case 'C' :
			case 'c' :
				return allecount[1];
			case 'G' :
			case 'g' :
				return allecount[2];
			case 'T' :
			case 't' :
				return allecount[3];
			default :
				fprintf(stderr,"Unrecognized allele %c in call to SiteData::cov\n",allele);
				return 0;
		}
	}
	else
		return depth;
}

std::string SiteData::id() const
{
	if (!_id.empty())
		return _id;
	else
		return "";
}

std::string& SiteData::id()
{
	return _id;
}

// Pileup class constructor
Pileup::Pileup ()
	: _fail(0),
	  _name(""),
	  _pos(0),
	  _refallele('\0'),
	  _encode(33.0),
	  _minQ(13.0),
	  _nind(0),
	  _depthReserve(20),
	  _numalt(0),
	  _numref(0),
	  _alleles(),
	  _poolsz(2)
{
	for (int i=0; i<2; ++i)
		_majmin[i]='\0';
}

// Pileup::setQualCode assigns quality encoding
void Pileup::setQualCode (double code)
{
	if (code < 0.0)
	{
		fprintf(stderr, "Invalid quality score encoding '%f' in Pileup::getCode\n", code);
		_fail = 1;
	}
	_encode = code;
}

// Pileup::setMinQ assigns minimum quality score for reads to be considered
void Pileup::setMinQ (double q)
{
	if (q < 0.0)
	{
		fprintf(stderr, "Cannot set minimum quality score to < 0.0 in Pileup::setMinQ\n");
		_fail = 1;
	}
	_minQ = q;
}

double Pileup::getMinQ ()
{
 return _minQ;
}

double Pileup::getQualCode ()
{
	return _encode;
}

// Pileup::getSeqDat extracts information from pileup line
int Pileup::getSeqDat (const std::string& pile, const char delim)
{
	int i;
	std::vector<std::string> pvec;
	tokenizeLine(pile, &pvec, delim);

	// extract sequence id, position, and reference allele
	setBasicSiteInfo(pvec);

	// initialize sequencing data vector if needed
	if (seqdat.empty())
	{
		_nind = extractIndN(pile, delim);
		if (_nind < 1)
		{
			_fail = 1;
			throw PileupException((ExceptionFormatter() << "No individuals found in pileup input in call to Pileup::" << __func__ << "()").str().c_str());
		}
		initializeSeqdat(_nind);
	}

	// initialize values
	_numref = 0;
	_numalt = 0;
	for (i = 0; i < 5; ++i)
	{
		_alleles[i] = 0;
		_wtalleles[i] = 0.0;
	}
	for (i=0; i<2; ++i)
		_majmin[i] = '\0';
	for(unsigned int k=0; k < _tcounts.size(); ++k)
	{
		_tcounts[k].total = 0.0;
		for (int j=0; j<5; ++j)
			_tcounts[k].counts[j] = 0.0;
	}

	// set reads and quality scores
	unsigned int ind = 0;
	std::string reads;
	std::string qscores;
	std::vector<std::string>::const_iterator iter = pvec.begin() + 3;

	if ((*iter).empty())
	{
		_fail = 1;
		throw PileupFormatException((ExceptionFormatter() << "No sequencing data for " << _name << " " << _pos << " in call to Pileup::" << __func__ << "()").str().c_str());
	}

	for (iter = pvec.begin() + 3; iter != pvec.end(); ++iter)
	{
		seqdat[ind].depth = 0;
		for (int i = 0; i < 4; ++i)
			seqdat[ind].allecount[i] = 0;
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
			getReadDat(&reads, &qscores, ind);
		}
		++ind;
	}
	return 0;
}

void Pileup::setBasicSiteInfo (const std::vector<std::string>& linevec)
{
	// set name
	if (!linevec[0].empty())
		_name = linevec[0];
	else
	{
		_fail = 1;
		throw PileupFormatException("Attempt to assign empty name in Pileup::setBasicSiteInfo()");
	}

	// set position
	if (!linevec[1].empty())
		_pos = atoi(linevec[1].c_str());
	else
	{
		_fail = 1;
		throw PileupFormatException("Attempt to assign empty position in Pileup::setBasicSiteInfo()");
	}

	// set reference allele identity
	if (!linevec[2].empty())
		_refallele = static_cast<char>(toupper(linevec[2][0]));
	else
	{
		_fail = 1;
		throw PileupFormatException("Attempt to assign empty reference allele in Pileup::setBasicSiteInfo()");
	}
}

void Pileup::tokenizeLine (const std::string& s, std::vector<std::string>* elems, const char delim)
{
	std::stringstream ss(s);
	std::string sholder;
	while (std::getline(ss, sholder, delim))
	{
		if (!sholder.empty())
			elems->push_back(sholder);
	}

	if (elems->size() < 3)
		throw PileupFormatException("Pileup::tokenizeLine() unable to split line");
}

// assigns read and quality scores for an individual
void Pileup::getReadDat (std::string* reads, std::string* qual, unsigned int ind)
{
	enum readtype_t {A_READ, C_READ, G_READ, T_READ, INDEL_READ, N_READ};
	readtype_t refidx = N_READ;
	char r;
	unsigned int index = 0;
	unsigned int indsize = 0;
	unsigned int read_num = 0;
	std::string::iterator q = qual->begin();

	switch (_refallele)
	{
                case 'A' :
                        refidx=A_READ;
                        break;
                case 'C' :
                        refidx=C_READ;
                        break;
                case 'G' :
                        refidx=G_READ;
                        break;
                case 'T' :
                        refidx=T_READ;
                        break;
                case 'N' :
                		refidx=N_READ;
                		break;
                default :
                	fprintf(stderr, "Warning: unrecognized reference allele '%c' at %s position %d\n", _refallele, _name.c_str(), _pos);
	}

	for(std::string::const_iterator r_iter = reads->begin(); r_iter != reads->end(); ++r_iter) // go over all reads
	{
		r = static_cast<char>(toupper(*r_iter));

		switch(r)
		{
			/* reference allele */
			case '.' :
			case ',' :
				try
				{
					recordRead(ind, q, read_num, index, refidx, _refallele);
				}
				catch (UnknownReadException& error)
				{
					std::cerr << error.what() << ": " << _name << " " << _pos << ", individual " << ind << " -> skipping read\n";
					++q;
					++index;
				}
				break;
			/* alternate allele */
			case 'A' :
				recordRead(ind, q, read_num, index, A_READ, r);
				break;
			case 'C' :
				recordRead(ind, q, read_num, index, C_READ, r);
				break;
			case 'G' :
				recordRead(ind, q, read_num, index, G_READ, r);
				break;
			case 'T' :
				recordRead(ind, q, read_num, index, T_READ, r);
				break;
			/* missing read or gap */
			case 'N' :
			case '*' :
			case '<' :
			case '>' :
				++q;
				++index;
				break;
			/* indel */
			case '+' :
			case '-' :
				indsize = indelSize(reads, (index+1));
                r_iter += indsize;
                index += indsize + 1;
                ++_alleles[INDEL_READ];
				break;
			/* beginning of read */
			case '^' :
				++r_iter;
				index += 2;
				break;
			/* end of read */
			case '$' :
				++index;
				break;
			/* unrecognized symbol */
			default :
				fprintf(stderr, "Unrecognized symbol '%c' in %s position %u:\n%s\t%s",r,_name.c_str(),_pos,reads->c_str(),qual->c_str());
				_fail = 1;
		}

	}
	seqdat[ind].depth = read_num;
}

void Pileup::recordRead (const unsigned int ind, std::string::iterator& q, unsigned int& read_num, unsigned int& index, const int id, const char read)
{
	// check for invalid read id
	if (id > 4) throw UnknownReadException(read); // A=>0, C=>1, G=>2, T=>3, INDEL=>4

	const float factor = 0.5;
	bool add = true;
	double phredq = static_cast<double>(*q) - _encode;
	double wtcount;

	if (phredq < 0.0)
	{
		_fail = 1;
		throw std::runtime_error(ExceptionFormatter() << "Negative quality score " << phredq << " at " << _name << " " << _pos << " in call to Pileup::" << __func__ << "()");
	}

	if (phredq >= _minQ)
	{
        	if (read_num >= seqdat[ind].rdat.size())
        	{
                	if (read_num < seqdat[ind].rdat.capacity())
                	{
                        	seqdat[ind].rdat.push_back(std::make_pair(read, phredq));
                        	add=false;
                	}
                	else
                        	seqdat[ind].rdat.resize(seqdat[ind].rdat.size() + static_cast<size_t>(factor * static_cast<float>(seqdat[ind].rdat.size())));
		}
		if (add)
		{
        		seqdat[ind].rdat[read_num].first = read;
        		seqdat[ind].rdat[read_num].second = phredq;
		}
		++read_num;
		read == _refallele ? ++_numref : ++_numalt;
		++seqdat[ind].allecount[id];
		++_alleles[id];
		wtcount = 1.0 - error(phredq);
		_wtalleles[id] += wtcount;
		addtreatcount(wtcount, &seqdat[ind]._id, id);
	}

	++q;
	++index;
}

void Pileup::addtreatcount (double amount, std::string* id, int allele)
{
	/*
	 * allele: 0=>A, 1=>C, 2=>G, 3=>T, 4=>INDEL
	 */

	static std::map<std::string, int> indexmap;
	unsigned int i;
	if (_tcounts.size() > 0 && id->size() > 0)
	{
		if (indexmap.find(*id) != indexmap.end())
		{
			i = indexmap[*id];
			_tcounts[i].counts[allele] += amount;
			_tcounts[i].total += amount;
		}
		else
		{
			for (i=0; i < _tcounts.size(); ++i)
			{
				if (_tcounts[i].name == *id)
				{
					indexmap[*id] = i;
					break;
				}
			}
			if (indexmap.find(*id) != indexmap.end())
			{
				_tcounts[i].counts[allele] += amount;
				_tcounts[i].total += amount;
			}
		}
	}
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
        return static_cast<unsigned int>((atoi(indsize.c_str()) + indsize.length()));
}

// Pileup::fail returns value of _fail member
int Pileup::fail ()
{
	return _fail;
}

unsigned int Pileup::extractIndN (const std::string& line, const char delim)
{
	unsigned int nind = 0;
	int inc = 0;
	std::vector<std::string> ptoke;
	tokenizeLine(line, &ptoke, delim);

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
				++inc;
            iter += 2 + inc;
		}
	}
	return nind;
}

void Pileup::setIndN (unsigned int n)
{
	if (n > 0)
		_nind = n;
	else
		fprintf(stderr, "Attempt to set nonpositive number of individuals in Pileup::%s()\n", __func__);
}

void Pileup::initializeSeqdat (const size_t n)
{
	seqdat.resize(n);
	for (size_t i = 0; i < n; ++i)
	{
		seqdat[i].rdat.resize(_depthReserve);
		seqdat[i].depth = _depthReserve;
	}
	_nind=static_cast<unsigned int>(seqdat.size());
}

unsigned int Pileup::setn (const std::string& ins, const char delim)
{
	if (!ins.empty())
		initializeSeqdat(extractIndN(ins, delim));
	else
		throw PreConditionException((ExceptionFormatter() << "Empty pileup line passed to Pileup::" << __func__ << "()").str().c_str());

	return static_cast<unsigned int>(seqdat.size());
}

std::string Pileup::seqName () const
{
	return _name;
}

unsigned int Pileup::position () const
{
	return _pos;
}

void Pileup::setDepthReserve (unsigned int depth)
{
	_depthReserve = depth;
}

size_t Pileup::nInd () const
{
	return _nind;
}

size_t Pileup::altcount () const
{
	return _numalt;
}

size_t Pileup::refcount () const
{
	return _numref;
}

double Pileup::altfreq () const
{
	return  static_cast<double>(_numalt) / (_numref + _numalt);
}

size_t Pileup::siteDepth () const
{
	return _numref + _numalt;
}

unsigned int Pileup::alleleCount (const char allele) const
{
	unsigned int n = 0;
	switch (toupper(allele))
	{
		case 'A' :
			n = _alleles[0];
			break;
		case 'C' :
			n = _alleles[1];
			break;
		case 'G' :
			n = _alleles[2];
			break;
		case 'T' :
			n = _alleles[3];
			break;
		case 'I' :
			n = _alleles[4];
			break;
		default :
			fprintf(stderr, "Unrecognized base '%c' in call to Pileup::alleleCount\n", allele);
	}
	return n;
}

double Pileup::wtalleleCount (const char allele) const
{
	double n = 0;
	switch (toupper(allele))
	{
		case 'A' :
			n = _wtalleles[0];
			break;
		case 'C' :
			n = _wtalleles[1];
			break;
		case 'G' :
			n = _wtalleles[2];
			break;
		case 'T' :
			n = _wtalleles[3];
			break;
		case 'I' :
			n = _wtalleles[4];
			break;
		default :
			fprintf(stderr, "Unrecognized base '%c' in call to Pileup::wtalleleCount\n", allele);
	}
	return n;
}

char Pileup::refAllele () const
{
	return _refallele;
}

void Pileup::assignTreatment (std::string id)
{
	_treatment = id;
}

std::string Pileup::treatment () const
{
	return _treatment;
}

void Pileup::setpoolsz (unsigned int n)
{
	/* n is the haploid sample size of each pool
	 * for diploids, n=2
	 */
	_poolsz = n;
}

unsigned int Pileup::poolsz (int ploidy) const
{
	/*
	 * if ploidy argument is supplied, number of individuals comprising each pool returned
	 * otherwise, the haploid size of each pool is returned
	 */

	if (_poolsz % ploidy == 0)
		return _poolsz/ploidy;
	else
		fprintf(stderr,"Incorrect ploidy passed to Pileup::poolsz(): haploid size = %u, ploidy = %i\n",_poolsz,ploidy);
	return 0;
}

void Pileup::setMajor(char allele)
{
	allele = static_cast<char>(toupper(allele));
	switch (allele)
	{
		case 'A' :
		case 'C' :
		case 'G' :
		case 'T' :
		case 'N' :
			_majmin[0] = allele;
			break;
	default :
		fprintf(stderr,"Invalid allele '%c' passed to Pileup::%s()\n",allele,__func__);
	}
}

void Pileup::setMinor(char allele)
{
	allele = static_cast<char>(toupper(allele));
	switch (allele)
	{
		case 'A' :
		case 'C' :
		case 'G' :
		case 'T' :
		case 'N' :
			_majmin[1] = allele;
			break;
	default :
		fprintf(stderr,"Invalid allele '%c' passed to Pileup::%s()\n",allele,__func__);
	}
}

char Pileup::majorid () const
{
	return _majmin[0];
}

char Pileup::minorid () const
{
	return _majmin[1];
}

char Pileup::empiricalMajor (bool wt)
{
	char maj='A';
	static const char a[] = {'A', 'C', 'G', 'T'};
	for (int i=0; i<4; ++i)
	{
		if (!wt)
		{
			if (_alleles[i] > alleleCount(maj))
				maj=a[i];
		}
		else
		{
			if (_wtalleles[i] > wtalleleCount(maj))
				maj=a[i];
		}
	}
	return maj;
}

char Pileup::empiricalMinorFast (bool wt)
{
	int i=0;
	static const char a[] = {'A', 'C', 'G', 'T'};
	char minor = a[i];
	char maj = empiricalMajor(wt);
	while (a[i] == maj && i < 4) minor = a[++i];
	for (i=0; i<4; ++i)
	{
		if (!wt)
		{
			if (_alleles[i] > alleleCount(minor) && a[i] != maj) minor = a[i];
		}
		else
		{
			if (_wtalleles[i] > wtalleleCount(minor) && a[i] != maj) minor = a[i];
		}
	}
	return minor;
}

bool Pileup::countCmp (std::pair<char,double> i, std::pair<char,double> j)
{
	return (i.second < j.second);
}


bool Pileup::baseCmp (std::pair<char,double> i, std::pair<char,double> j)
{
	return ((int)i.first < (int)j.first);
}

char Pileup::empiricalMinor ()
{
	static const seqread b [] = {std::make_pair('A',0.0), std::make_pair('C',0.0), std::make_pair('G',0.0), std::make_pair('T',0.0)};
	static std::vector<seqread> counts (b, b + sizeof(b)/sizeof(seqread));
	static std::vector<seqread>::iterator i;
	static std::vector<SiteData>::iterator j;

	// sort the elements of counts vector in A,C,G,T order and set counts to zero
	std::sort (counts.begin(), counts.end(), baseCmp);

	for (i=counts.begin(); i!=counts.end(); ++i)
		i->second=0.0;

	// count occurrence of each base at site weighted by the quality score
	for (j=seqdat.begin(); j!=seqdat.end(); ++j)
	{
		for (unsigned int k=0; k<j->cov(); ++k)
		{
			switch (j->rdat[k].first)
			{
				case 'A' :
					counts[0].second += 1.0-error(j->rdat[k].second);
					break;
				case 'C' :
					counts[1].second += 1.0-error(j->rdat[k].second);
					break;
				case 'G' :
					counts[2].second += 1.0-error(j->rdat[k].second);
					break;
				case 'T' :
					counts[3].second += 1.0-error(j->rdat[k].second);
					break;
			}
		}
	}

	// sort the elements of counts vector in ascending count order and find second most common base
	std::sort (counts.begin(), counts.end(), countCmp);
	if (counts[2].second == 0.0)
		return 'N'; // site is fixed, i.e. no minor allele
	else
		return counts[2].first;
}

double Pileup::error (double q)
{
	const int n = 71;
	static double err[n];
	static int qidx = static_cast<int>(q);

	if (q < 0.0)
	{
		fprintf(stderr,"Invalid quality score '%f' passed to Pileup:error\n",q);
		_fail = 1;
		return 1.0;
	}

	if (qidx < n)
		return (err[qidx] != 0.0 ? err[qidx] : (err[qidx]=scalePhred(q)));
	else
		return scalePhred(q);
}

double Pileup::scalePhred (double q)
{
	return pow(10, -q/10.0);
}

void Pileup::addtreatment (std::string* name)
{
	_tcounts.resize(_tcounts.size() + 1);
	for (unsigned int i=0; i<_tcounts.size(); ++i)
	{
		if (_tcounts[i].name.size() == 0)
		{
			_tcounts[i].name = *name;
			for (int j=0; j<5; ++j)
				_tcounts[i].counts[j] = 0.0;
			break;
		}
	}
}

double Pileup::treatcounts (std::string* id, char allele)
{
	static std::map<std::string, int> index;
	unsigned int i;
	double c = 0.0;
	if (index.find(*id) != index.end())
	{
		i = index[*id];
		if (!allele)
			return _tcounts[i].total;
		c = tcount(i, allele);
	}
	else
	{
		for(i=0; i < _tcounts.size(); ++i)
		{
			if (_tcounts[i].name == *id)
			{
				index[*id] = i;
				break;
			}
		}
		if (index.find(*id) != index.end())
		{
			if (!allele)
				return _tcounts[i].total;
			c = tcount(i, allele);
		}
	}
	return c;
}

double Pileup::tcount (int i, char a)
{
	switch (a)
	{
		case 'A' :
			return _tcounts[i].counts[0];
		case 'C' :
			return _tcounts[i].counts[1];
		case 'G' :
			return _tcounts[i].counts[2];
		case 'T' :
			return _tcounts[i].counts[3];
		case 'I' :
			return _tcounts[i].counts[4];
	}
	return 0.0;
}

unsigned int Pileup::idSize (const std::string& line, const char delim)
{
	std::vector<std::string> tokens;
	tokenizeLine(line, &tokens, delim);
	return static_cast<unsigned int>(tokens[0].size());
}

PileupException::PileupException(const char* error)
	: error_(error),
	  _errmsg("")
{}

PileupException::~PileupException () throw() {}

const char* PileupException::what() throw()
{
	std::stringstream msg;
	msg << "Pileup exception occurred:\n";
	if (error_) msg << error_;
	if (!_errmsg.empty()) _errmsg.clear();
	_errmsg = msg.str();
	return _errmsg.c_str();
}

PileupFormatException::PileupFormatException(const char* error)
	: PileupException(error)
{}

const char* PileupFormatException::what() throw()
{
	std::stringstream msg;
	msg << PileupException::what() << "\nInvalid pileup file format\n";
	if (!_errmsg.empty()) _errmsg.clear();
	_errmsg = msg.str();
	return _errmsg.c_str();
}

UnknownReadException::UnknownReadException (const char readtype)
	: std::runtime_error(ExceptionFormatter() << readtype << " is an invalid read type") {}
