//parsePileup.h

#ifndef PARSEPILEUP_H_
#define PARSEPILEUP_H_

#include <vector>
#include <string>
#include <utility>

// STRUCTURE DEFINITIONS
struct SiteData
{
public:
	SiteData ();
	std::vector< std::pair<int, int> > rdat;
	size_t depth;
private:
};

// CLASS DEFINITIONS
class Pileup
{
public:
	Pileup ();
	std::string seqName ();
	unsigned int position ();
	int getSeqDat (const std::string&); // extracts data from pileup line
	void setQualCode (int); // sets encode
	void setMinQ (int q); // sets minimum quality score for read to be considered
	int getMinQ (); // return minQ value
	int getQualCode (); // returns encode value
	size_t ExtractIndN (const std::string& line); // find number of individuals in pileup line
	void setIndN (size_t n); // sets _nind member
	void setDepthReserve (size_t depth);
	size_t altcount (); // return number of alternate alleles for site
	size_t refcount (); // return number of reference alleles for site
	double altfreq (); // return alternate allele frequency for site
	size_t siteDepth (); // return total number reads for site
	size_t nInd ();
	unsigned int alleleCount (char allele); // return count of specific alleles: 'A', 'C', 'G', 'T', 'I' (indel)
	char refAllele (); // return reference allele
	int fail (); // returns value of _fail member
	std::vector<SiteData> seqdat; // stores reads and quality scores
private:
	void getReadDat (std::string*, std::string*, size_t);
	unsigned int indelSize (std::string*, unsigned int);
	void initializeSeqdat (size_t n); // intializes seqdat member
	void newRead (const size_t ind, const size_t read_num, int read, int quality);
	std::string name; // chromosome name
	unsigned int pos; // position
	int encode; // quality score encoding: 33 for illumina v1.8+ or 64 for illumina v1.3+
	int minQ; // minimum quality score for read to be kept
	int _fail;
	size_t _nind;
	size_t _depthReserve;
	size_t _numalt;
	size_t _numref;
	unsigned int _alleles[5]; // counts of specific alleles for site [A,C,G,T,indel]
	char _refallele; // reference allele
};

#endif /* PARSEPILEUP_H_ */
