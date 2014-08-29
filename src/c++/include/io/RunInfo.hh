/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Component to read/write RunInfo.xml files.
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_IO_RUNINFO_HH
#define EAGLE_IO_RUNINFO_HH

#include <vector>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/utility.hpp>

using namespace std;


namespace eagle
{
namespace io
{

struct ReadDescription
{
    ReadDescription() : firstCycle(0), lastCycle(0), isIndex(false) {}
    /*
public:
    ReadDescription( const ReadDescription &obj) : firstCycle(obj.firstCycle), lastCycle(obj.lastCycle), isIndex(obj.isIndex) {}
    */
    unsigned int firstCycle, lastCycle;
    bool isIndex;
};

class RunInfo
{
public:
    RunInfo( const boost::filesystem::path &filename ) 
        : runId("")
        , runNumber("")
        , tileNameMethod("")
        , flowcell("")
        , laneCount(0)
        , surfaceCount(0)
        , swathCount(0)
        , tileCount(0)
        , reads()
    { 
        parse(filename);
    }
    void parse( const boost::filesystem::path &filename );

    unsigned int getClusterLength() const
    {
        assert (!reads.empty());
        return reads.rbegin()->lastCycle;
    }

    string runId;
    string runNumber;
    string tileNameMethod;
    string flowcell;
    unsigned int laneCount, surfaceCount, swathCount, tileCount;
    std::vector<ReadDescription> reads;
};


} // namespace io
} // namespace eagle

#endif // EAGLE_IO_RUNINFO_HH
