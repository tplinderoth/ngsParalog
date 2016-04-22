//parsePileup.h

#ifndef PARSEPILEUP_H_
#define PARSEPILEUP_H_

#include <vector>
#include <string>
#include <utility>

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
	size_t depth;
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
class Pileup
{
public:
	Pileup ();
	std::string seqName () const;
	unsigned int position () const;
	int getSeqDat (const std::string&); /* extracts data from pileup line */
	void setQualCode (double); /* sets encode, the phred score offset value */
	void setMinQ (double q); /* sets minimum quality score for read to be considered */
	int getMinQ (); /* return minQ value */
	int getQualCode (); /* returns encode value */
	size_t ExtractIndN (const std::string& line); /* find number of individuals in pileup line */
	void setIndN (size_t n); /* sets _nind member */
	void setDepthReserve (size_t depth);
	size_t altcount () const; /* return number of alternate alleles for site */
	size_t refcount () const; /* return number of reference alleles for site */
	double altfreq () const; /* return alternate allele frequency for site */
	size_t siteDepth () const; /* return total number reads for site */
	size_t nInd () const;
	void setpoolsz (unsigned int n); /* sets _poolsz member */
	unsigned int poolsz (int ploidy=1) const; /* if ploidy argument is supplied, number of individuals comprising each pool returned, otherwise return haploid size*/
	unsigned int alleleCount (const char allele) const; /* return count of specific alleles: 'A', 'C', 'G', 'T', 'I' (indel) for site*/
	double wtalleleCount (const char allele) const; /* return quality-score-weighted count of a specific allele for site*/
	char refAllele () const; /* return reference allele */
	void assignTreatment (std::string id); /* sets _treatment member */
	std::string treatment () const; /* returns treatment member */
	unsigned int setn (std::string ins); /* extracts number of individuals from pileup line and initializes seqdat member*/
	void setMajor(char allele); /* assign major allele*/
	void setMinor (char allele); /* assign minor allele */
	char majorid () const; /* returns major allele. If not previously set, returns empirical major */
	char minorid () const; /* returns minor allele. If not previously set, returns empirical minor */
	char empiricalMajor (bool wt=false); /* returns the most common base for site */
	char empiricalMinor (); /* returns second most common base - counts are weighted by quality scores */
	char empiricalMinorFast (bool wt=false); /* returns second most common base without considering quality scores */
	double error (int q); /* returns probability of error indicated by base quality score */
	int fail (); /* returns value of _fail member */
	static bool countCmp (seqread i, seqread j); /* compares the double values of seqread type */
	static bool baseCmp (seqread i, seqread j); /* compares the ascii encoded decimal values of the char value for seqread type */
	double treatcounts (std::string* id, char allele = '\0'); /* returns count of allele in treatment 'id' or total count if allele not specified */
	void addtreatcount (double amount, std::string* id, int allele); /* add counts to _tcounts 0=A,1=C,2=G,3=T,4=INDEL*/
	void addtreatment (std::string* name); /* adds a TreatCounts element to _tcounts */
	std::vector<SiteData> seqdat; /* stores reads and quality scores */
private:
	void getReadDat (std::string*, std::string*, size_t);
	unsigned int indelSize (std::string*, unsigned int);
	void initializeSeqdat (size_t n); /* intializes seqdat member */
	void recordRead (const size_t ind, std::string::iterator& q, size_t& read_num, unsigned int& index, const int id, const char read);
	double scalePhred (double q); /* rescales phred quality score */
	double tcount (int i, char a); /* used with treatcounts to get counts from _tcounts */
	std::string _name; /* chromosome name */
	unsigned int _pos; /* position */
	double _encode; /* minimum possible ASCII decimal value used to encode quality scores */
	double _minQ; /* minimum quality score for read to be kept */
	int _fail;
	size_t _nind;
	size_t _depthReserve;
	size_t _numalt;
	size_t _numref;
	double _alleles[5]; /* counts of specific alleles for site [A,C,G,T,INDEL] */
	double _wtalleles[5]; /* counts of specific alleles weighed by quality scores [A,C,G,T,INDEL]*/
	char _refallele; /* reference allele */
	std::string _treatment; /* associates pileup object with a particular treatment identifier (subset of individuals) */
	unsigned int _poolsz; /* haploid sample size constituting each pool (i.e. element of seqdat member) */
	char _majmin [2]; /* [major allele, minor allele] */
	std::vector<TreatCounts> _tcounts; /* counts of alleles for each treatment weighted by the quality scores */
};

#endif /* PARSEPILEUP_H_ */
