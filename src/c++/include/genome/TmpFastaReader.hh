/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Fast fasta reader, which we keep until the mainstream reader is faster than this one
 **
 ** \author Lilian Janin
 **/


#ifndef TMP_FASTA_READER_HH
#define TMP_FASTA_READER_HH

#include <fstream>
#include <vector>
#include <boost/filesystem.hpp>

#include "genome/Reference.hh"

using namespace std;


namespace eagle
{
namespace genome
{


struct TmpFastaFileInfo
{
    ifstream *file;
    unsigned long globalPosMin, globalPosMax;
    unsigned int headerLength;
    unsigned int basesPerLine;
    unsigned long baseCount;
    string contigName;
    unsigned long fileSize;
};


class TmpFastaReader
{
public:
    TmpFastaReader( const boost::filesystem::path &refDir );
    const vector<string> allContigNames() { return fastaRef_.allContigNames(); }
    const vector<unsigned long> allContigLengths() { return fastaRef_.allContigLengths(); }
    char get( const unsigned long globalPos, const unsigned long offset, bool& overlapContigBoundary );
    void convertFromGlobalPos( const unsigned long globalPos, int& refId, unsigned long& posInContig );

private:
    MultiFastaReference fastaRef_;
    vector<TmpFastaFileInfo> fileInfos_;
};


} // namespace genome
} // namespace eagle

#endif
