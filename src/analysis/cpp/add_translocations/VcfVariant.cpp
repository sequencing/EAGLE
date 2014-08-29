#include "VcfVariant.h"


void VcfParser::buildVcfRecordFromString(const std::string& vcfLine, VcfRecord& rec) {
	std::vector<std::string> vcfTokens;
	_tokenizeLine(vcfLine,vcfTokens);
	switch (vcfTokens.size()) {
		case 10:
			rec.chr=vcfTokens[0]; rec.pos=_stringToInt(vcfTokens[1]); rec.id=vcfTokens[2]; rec.ref=vcfTokens[3];
			rec.alt=vcfTokens[4]; rec.qual=_stringToInt(vcfTokens[5]); rec.filter=vcfTokens[6]; rec.info=vcfTokens[7];
			rec.format=vcfTokens[8]; rec.gt=vcfTokens[9];
			break;
		default:
			std::cout << "Cannot parse line " << vcfLine << std::endl;
			break;
	}
	bool reverse(false);
	size_t idx1 = rec.id.find("bnd_U");
	size_t idx2 = rec.id.find("bnd_X");
	if (idx1 != std::string::npos || idx2 != std::string::npos) {
		reverse=true;
	}

	if (!_parseInfoFieldAndAssignLength(rec,reverse)) {
		std::cerr << "WARNING : cannot parse variant length from variant info field : " << rec.info << std::endl;
	}
}

bool VcfParser::_parseInfoFieldAndAssignLength(VcfRecord& rec, const bool reverse) {
	bool svLenFound(false);
	boost::tokenizer< boost::char_separator<char> > infoTokens(rec.info, boost::char_separator<char>(";"));
	BOOST_FOREACH (const std::string& t, infoTokens) {
		//std::cout << "token=" << t << std::endl;
		size_t idx = t.find("SVLEN");
		if (idx != std::string::npos) {
			//std::cout << "SVLEN found in: " << t << std::endl;
			svLenFound = true;
			std::vector<std::string > fields; 
			boost::split(fields,t,boost::is_any_of("="));
			try {
				rec.len = abs(boost::lexical_cast<int>(fields[1])); 
				if (reverse) {
					rec.len *= -1;
				}
			} catch (boost::bad_lexical_cast &) {
				std::cout << "Error while trying to convert SVLEN info entry to numerical value: " << t << std::endl;
			}
			//std::cout << "Extracted event length=" << rec.len << std::endl;
			//assert(rec.len>0);
			svLenFound = true;
		}
	}
	return svLenFound;
}

void VcfParser::_tokenizeLine(const std::string& vcfLine, std::vector<std::string>& vcfTokens) {
	boost::tokenizer< boost::char_separator<char> > vcfLineTokens(vcfLine,boost::char_separator<char>("\t"));
	BOOST_FOREACH(const std::string& t, vcfLineTokens) { vcfTokens.push_back(t); }
}

int VcfParser::_stringToInt( const std::string& s ) {
	std::istringstream i(s);
	int x;
	if (!(i >> x))
		return 0;
	return x;
} 

std::ostream & operator<<(std::ostream& os, const VcfRecord & r) {
	return os << r.chr << "\t" << r.pos << "\t" << r.id << "\t" << r.ref << "\t" 
            << r.alt << "\t" << r.qual << "\t" << r.filter << "\t" << r.info 
            << "\t" << r.format << "\t" << r.gt;
}

