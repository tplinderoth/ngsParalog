//parsePileup.h

#ifndef PARSEPILEUP_H_
#define PARSEPILEUP_H_

#include <vector>
#include <string>
#include <utility>
#include <stdexcept>

typedef std::pair<char,double> seqread;

// STRUCTURE DEFINITIONS

struct SiteData
{
friend class Pileup;
public:
	SiteData ();
	unsigned int cov (char allele='\0') const; /* returns count of 'allele' or total coverage if no argument is provided */
	std::string id() const; /* returns _id member */
	std::string& id(); /* returns &_id */
	std::vector<seqread> rdat;
private:
	unsigned int depth;
	int allecount [5]; /* counts of [A, C, G, T, INDEL] alleles */
	std::string _id; /* treatment affiliation */
};

struct TreatCounts
{
friend class Pileup;
private:
	std::string name; // treatment name
	double counts [5]; // counts of [A,C,G,T,INDEL] alleles
	double total; // total counts
};

// CLASS DEFINITIONS
class PileupException : public std::exception
{
public:
	explicit PileupException (const char* error = NULL);
	virtual const char* what () const throw();
protected:
	const char* error_;
};

class PileupFormatException : public PileupException
{
public:
	explicit PileupFormatException(const char* error = NULL);
	virtual const char* what () const throw();
};

class UnknownReadException : public std::runtime_error
{
public:
	UnknownReadException (const char readtype);
};

class Pileup
{
public:
	Pileup ();
	void setQualCode (double); /* sets encode, the phred score offset value */
	void setMinQ (double q); /* sets minimum quality score for read to be considered */
	unsigned int setn (const std::string& ins, const char delim='\t'); /* extracts number of individuals from pileup line and initializes seqdat member*/
	void setMajor(char allele); /* assign major allele*/
	void setMinor (char allele); /* assign minor allele */
	void setpoolsz (unsigned int n); /* sets _poolsz member */
	void assignTreatment (std::string id); /* sets _treatment member */
	void addtreatcount (double amount, std::string* id, int allele); /* add counts to _tcounts 0=A,1=C,2=G,3=T,4=INDEL*/
	void addtreatment (std::string* name); /* adds a TreatCounts element to _tcounts */
	void setDepthReserve (unsigned int depth);
	std::string seqName () const;
	unsigned int position () const;
	int getSeqDat (const std::string&, const char delim='\t'); /* extracts data from pileup line */
	double getMinQ (); /* return minQ value */
	double getQualCode (); /* returns encode value */
	static unsigned int extractIndN (const std::string& line, const char delim = '\t'); /* find number of individuals in pileup line */
	size_t altcount () const; /* return number of alternate alleles for site */
	size_t refcount () const; /* return number of reference alleles for site */
	double altfreq () const; /* return alternate allele frequency for site */
	size_t siteDepth () const; /* return total number reads for site */
	size_t nInd () const;
	unsigned int poolsz (int ploidy=1) const; /* if ploidy argument is supplied, number of individuals comprising each pool returned, otherwise return haploid size*/
	unsigned int alleleCount (const char allele) const; /* return count of specific alleles: 'A', 'C', 'G', 'T', 'I' (indel) for site*/
	double wtalleleCount (const char allele) const; /* return quality-score-weighted count of a specific allele for site*/
	char refAllele () const; /* return reference allele */
	std::string treatment () const; /* returns treatment member */
	char majorid () const; /* returns major allele. If not previously set, returns empirical major */
	char minorid () const; /* returns minor allele. If not previously set, returns empirical minor */
	char empiricalMajor (bool wt=false); /* returns the most common base for site */
	char empiricalMinor (); /* returns second most common base - counts are weighted by quality scores */
	char empiricalMinorFast (bool wt=false); /* returns second most common base without considering quality scores */
	double error (double q); /* returns probability of error indicated by base quality score */
	int fail (); /* returns value of _fail member */
	static bool countCmp (seqread i, seqread j); /* compares the double values of seqread type */
	static bool baseCmp (seqread i, seqread j); /* compares the ascii encoded decimal values of the char value for seqread type */
	double treatcounts (std::string* id, char allele = '\0'); /* returns count of allele in treatment 'id' or total count if allele not specified */
	static void tokenizeLine (const std::string& s, std::vector<std::string>* elems, const char delim='\t'); /* tokenizes pileup line using 'delim' as the delimiter */
	static unsigned int idSize (const std::string& line, const char delim='\t'); /* returns length of sequence identifier from pileup line */
	std::vector<SiteData> seqdat; /* stores reads and quality scores */
protected:
	void setBasicSiteInfo (const std::vector<std::string>& linevec); /* splits pileup line & sets _name, _pos, _refallele */
	void setIndN (unsigned int n); /* sets _nind member */
	int _fail;
	std::string _name; /* chromosome name */
	unsigned int _pos; /* position */
	char _refallele; /* reference allele */
private:
	void getReadDat (std::string*, std::string*, unsigned int);
	unsigned int indelSize (std::string*, unsigned int);
	void initializeSeqdat (const size_t n); /* initializes seqdat member */
	void recordRead (const unsigned int ind, std::string::iterator& q, unsigned int& read_num, unsigned int& index, const int id, const char read);
	double scalePhred (double q); /* rescales phred quality score */
	double tcount (int i, char a); /* used with treatcounts to get counts from _tcounts */
	double _encode; /* minimum possible ASCII decimal value used to encode quality scores */
	double _minQ; /* minimum quality score for read to be kept */
	unsigned int _nind;
	unsigned int _depthReserve;
	unsigned int _numalt;
	unsigned int _numref;
	unsigned int _alleles[5]; /* counts of specific alleles for site [A,C,G,T,INDEL] */
	double _wtalleles[5]; /* counts of specific alleles weighed by quality scores [A,C,G,T,INDEL] */
	std::string _treatment; /* associates pileup object with a particular treatment identifier (subset of individuals) */
	unsigned int _poolsz; /* haploid sample size constituting each pool (i.e. element of seqdat member) */
	char _majmin [2]; /* [major allele, minor allele] */
	std::vector<TreatCounts> _tcounts; /* counts of alleles for each treatment weighted by the quality scores */
};

#endif /* PARSEPILEUP_H_ */
