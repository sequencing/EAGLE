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
#include "common/Exceptions.hh"


namespace eagle
{
namespace io
{


class FastqTile
{
public:
    FastqTile( const unsigned long long expectedReadCount, const unsigned int clusterLength, const std::string &filenameTemplate, const std::string &statsFilenameTemplate, const std::string &filterFilename, const std::string &clocsFilename, const std::string &controlFilename, const bool verbose=true );
    void addClusterToRandomLocation( const char *bufCluster, const bool isPassingFilter = true );
    void flushToDisk() const;

private:
    void writeFastqFile( const unsigned int cycle ) const;
    void writeStatsFile( const unsigned int cycle ) const;
    void writeFilterFile() const;
    void writeClocsFile() const;
    void writeControlFile() const;

    unsigned long long expectedReadCount_;
    unsigned int clusterLength_;
    std::string filenameTemplate_;
    std::string statsFilenameTemplate_;
    std::string filterFilename_;
    std::string clocsFilename_;
    std::string controlFilename_;
    //    vector<bool> usedLocations;
    std::vector<std::vector<unsigned int> > stats_;
    char *ramTile_;
    std::vector<char> passFilter_;
};



} // namespace io
} // namespace eagle

#endif // EAGLE_IO_FASTQ_HH
