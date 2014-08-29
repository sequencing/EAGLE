/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Writer component for BAM files.
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_GENOME_BAM_HH
#define EAGLE_GENOME_BAM_HH
#include <iostream>
#include <ostream>

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>

#include "config.h"
//#include "common/Debug.hh"
#include "common/Exceptions.hh"
//#include "flowcell/TileMetadata.hh"

#define ASSERT_MSG(x,y) assert((x) && y)


#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "common/Exceptions.hh"
#include "genome/ReadCluster.hh"
#include "io/RunInfo.hh"
#include "genome/BamAdapters.hh"
#include "genome/SharedFastaReference.hh"


namespace eagle
{
namespace genome
{

class BamOrMetadataOutput
{
public:
    BamOrMetadataOutput( const boost::filesystem::path outDir, eagle::io::RunInfo &runInfo, PreferredFastaReader* fastaReference = eagle::genome::SharedFastaReference::get() );
    ~BamOrMetadataOutput();
    void init( const boost::filesystem::path outDir );
    void add( eagle::genome::ReadClusterWithErrors& readClusterWithErrors );
    void addRebased( eagle::genome::ReadClusterWithErrors& readClusterWithErrors, const signed long globalPosShift, const unsigned long firstPosToProcess = 0, const unsigned long lastPosToProcess = 0, const bool dropLastBase = false, const genome::RefToSampleSegment& cigarModifierHelper = genome::RefToSampleSegment() );

private:
    eagle::io::RunInfo &runInfo_;
    PreferredFastaReader* fastaReference_;
//    ofstream simout_;
    boost::iostreams::filtering_ostream bamStream_;
    boost::iostreams::filtering_ostream bgzfStream_;
    std::list<EagleBamAlignmentAdapter> reorderedAlignmentsWindow_;

    bool updateLhsCIGAR( const vector< unsigned int >& reorderedCIGAR, vector< unsigned int >& softClippedCIGAR, genome::ReadClusterWithErrors& readClusterWithErrors, const genome::RefToSampleSegment& cigarModifierHelper, const unsigned long firstPosToProcess, const signed long GlobalPos, unsigned int& FLAG, unsigned long& globalPosAfterSoftClipping );
    bool updateRhsCIGAR( const vector< unsigned int >& reorderedCIGAR, vector< unsigned int >& softClippedCIGAR, genome::ReadClusterWithErrors& readClusterWithErrors, const genome::RefToSampleSegment& cigarModifierHelper, const unsigned long firstPosToProcess, const signed long GlobalFirstPos, const signed long GlobalLastPos, unsigned int& FLAG, unsigned long& globalPosAfterSoftClipping );
    void softClipCIGAR( const vector< unsigned int > &CIGAR, unsigned int clippingLength, vector< unsigned int > &softClippedCIGAR );
    void flushReorderedAlignmentsUntilPos( unsigned long POS );
    void addToReorderedAlignments( EagleBamAlignmentAdapter &eagleBamAlignmentAdapter );
};


} // namespace genome
} // namespace eagle

#endif // EAGLE_IO_BAM_HH
