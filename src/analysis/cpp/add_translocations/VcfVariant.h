#ifndef VCF_VARIANT_H
#define VCF_VARIANT_H

#include <iostream>
#include <sstream>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

using std::string;

struct VcfRecord; // forward declaration because I like to put the typedefs at the beginning

typedef std::vector<VcfRecord>       vcfStore; 

struct VcfRecord {

	VcfRecord() : pos(0),len(0),qual(0) {}

	VcfRecord(std::string c, int p, int l, std::string i, std::string r, std::string a, 
					  int q, std::string fi, std::string in, std::string fo, std::string g)
            : chr(c), pos(p), len(l), id(i), ref(r), alt(a), qual(q),
              filter(fi), info(in), format(fo), gt(g) {}
	
	std::string chr;
	int pos;
	int len; // length can be negative for backward loops in a breakpoint
	std::string id;
	std::string ref;
	std::string alt;
	int qual;
	std::string filter;
	std::string info;
	std::string format;
	std::string gt;
};

class VcfParser {
	public:
		VcfParser() {};

		void buildVcfRecordFromString(const std::string& vcfLine, VcfRecord& rec); 

	private:
		bool _parseInfoFieldAndAssignLength(VcfRecord& rec, bool reverse);
		void _tokenizeLine(const std::string& vcfLine, std::vector<std::string>& vcfTokens);
		int _stringToInt( const std::string& s );
};

std::ostream & operator<<(std::ostream& os, const VcfRecord & r);

#endif // VCF_VARIANT_H
