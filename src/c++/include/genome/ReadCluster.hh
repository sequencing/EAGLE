/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_GENOME_READ_CLUSTER_HH
#define EAGLE_GENOME_READ_CLUSTER_HH

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>

#include "genome/Reference.hh"
#include "io/RunInfo.hh"
#include "model/Nucleotides.hh"
#include "genome/QualityModel.hh"
#include "genome/EnrichedFragment.hh"
#include "genome/ReferenceToSample.hh"


namespace eagle
{
namespace genome
{


class ReadClusterSharedData
{
public:
    ReadClusterSharedData( const unsigned int clusterLength, const eagle::io::RunInfo &runInfo, const boost::filesystem::path& sampleGenomeDir, const std::vector<boost::filesystem::path>& qualityTableFiles, const boost::filesystem::path& mismatchTableFile, const boost::filesystem::path& homopolymerIndelTableFile, const boost::filesystem::path& motifQualityDropTableFile, const boost::filesystem::path& qqTableFile, const unsigned int userRandomSeed, const std::vector< std::string >& errorModelOptions );

    unsigned int clusterLength_;
    const eagle::io::RunInfo &runInfo_;
    ErrorModel errorModel_;
    unsigned int userRandomSeed_;
    std::vector<FragmentStructure> multiplexedFragmentStructures;
};


class ReadClusterWithErrors
{
public:
    ReadClusterWithErrors(ReadClusterSharedData &sharedData, const EnrichedFragment& eFragment, boost::shared_ptr<boost::mt19937> randomGen);
    const char *getBclCluster( bool generateCigar = false, const bool dropLastBase = false );
    const std::vector<unsigned int>& getCigar( unsigned int readNum, const bool dropLastBase = false );
    unsigned int getUsedDnaLength( unsigned int readNum, const bool dropLastBase = false );
    const string getNucleotideOrQualitySequenceForRead( unsigned int readNum, bool getNucleotides, bool revComp, const bool dropLastBase = false );

private:
    ReadClusterSharedData &sharedData_;
    boost::shared_ptr< boost::mt19937 > randomGen_;

public:
    const EnrichedFragment eFragment_;
    //    unsigned int id_;

private:
    bool lazyEvaluationDone_;
    string buf_;
    std::vector< std::vector< unsigned int > > cigar_;
    std::vector< unsigned int > usedDnaLength_;
//    unsigned int read1PosInBuf_, read1Length_, read2PosInBuf_, read2Length_;
};


class ReadClusterFactory
{
public:
    ReadClusterFactory( const eagle::io::RunInfo &runInfo, const boost::filesystem::path& sampleGenomeDir, const std::vector<boost::filesystem::path>& qualityTableFiles, const boost::filesystem::path& mismatchTableFile, const boost::filesystem::path& homopolymerIndelTableFile, const boost::filesystem::path& motifQualityDropTableFile, const boost::filesystem::path& qqTableFile, const unsigned int userRandomSeed, const std::vector< std::string >& errorModelOptions );
    ReadClusterWithErrors getReadClusterWithErrors( const eagle::model::Fragment &f );

private:
    const eagle::io::RunInfo &runInfo_;
    ReadClusterSharedData sharedData_;
};


} // namespace genome
} // namespace eagle

#endif // EAGLE_GENOME_READ_CLUSTER_HH
