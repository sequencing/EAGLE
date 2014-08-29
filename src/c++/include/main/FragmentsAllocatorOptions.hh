/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Command line options for 'allocateFragments'
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_OPTIONS_FRAGMENTS_ALLOCATOR_HH
#define EAGLE_OPTIONS_FRAGMENTS_ALLOCATOR_HH

#include <string>
#include <map>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace eagle
{
namespace main
{

class FragmentsAllocatorOptions : public eagle::common::Options
{
public:
    FragmentsAllocatorOptions();
private:
    std::string usagePrefix() const {return std::string("Usage:\n")
                                          + std::string("       allocateFragments [parameters] [options]");}
    void postProcess(boost::program_options::variables_map &vm);

public:
    //    unsigned long readCount;
    boost::filesystem::path sampleGenomeDir;
    boost::filesystem::path outputDir;
    float coverageDepth;
    unsigned long tileCount;
    unsigned int basesPerCluster;
    std::string tls;
    struct {
        double min, median, max, lowStdDev, highStdDev;
        std::string M0, M1;
    } templateLengthStatistics;
    bool uniformCoverage;
    std::string tileAllocationMethodStr;
    enum { TILE_ALLOCATION_RANDOM, TILE_ALLOCATION_SEQUENCE, TILE_ALLOCATION_INTERLEAVED } tileAllocationMethod;
    unsigned int randomSeed;
    boost::filesystem::path templateLengthTableFile;
    std::string contigName;
    bool mergeExistingFragments;
    boost::filesystem::path gcCoverageFitFile;
    double maxCoverageError;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_OPTIONS_FRAGMENTS_ALLOCATOR_HH
