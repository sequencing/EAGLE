/**
 ** Copyright (c) 2018 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Writer component for FASTQ files.
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_IO_FASTQ_HH
#define EAGLE_IO_FASTQ_HH


#include <string>
#include <vector>
#include <iostream>
#include "common/Exceptions.hh"
#include "io/RunInfo.hh"


namespace eagle
{
namespace io
{


class FastqTile
{
public:
    FastqTile( const unsigned long long expectedReadCount, const unsigned int clusterLength, const std::string &read1FastqFilename, const std::string &read2FastqFilename, const RunInfo &runInfo, const int lane, const unsigned int tileId, const bool verbose=true );

    void addCluster( const std::string &read1Nucleotides, const std::string &read1Qualities, const std::string &read2Nucleotides, const std::string &read2Qualities, const bool isPassingFilter = true );
    void finaliseAndWriteInfo();

private:

    unsigned long long expectedReadCount_;
    unsigned int clusterLength_;
    std::string filenameTemplate_;
    std::string read1FastqFilename_;
    std::string read2FastqFilename_;
    std::ofstream read1FastqFile_;
    std::ofstream read2FastqFile_;
    std::ofstream infoFile_;

    std::string readNamePrefix_;
    unsigned long long totalReadCount_;
    unsigned long long passedFilterReadCount_;
};



} // namespace io
} // namespace eagle

#endif // EAGLE_IO_FASTQ_HH
