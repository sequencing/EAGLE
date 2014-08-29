/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MAIN_BAM_ANALYSER_HH
#define EAGLE_MAIN_BAM_ANALYSER_HH

#include <vector>
#include <boost/iostreams/filtering_stream.hpp>

#include "genome/Reference.hh"
#include "io/BamParserFilter.hh"
#include "io/StorableBamAlignment.hh"
#include "model/FragmentPosResolver.hh"
#include "BamDiffOptions.hh"


namespace eagle
{
namespace main
{


class BamDiff
{
public:
    BamDiff (const BamDiffOptions &options);
    ~BamDiff();
    void run();

private:
    void CreateBamOutputStream( const boost::filesystem::path& outFilename, boost::iostreams::filtering_ostream& bamStream_, boost::iostreams::filtering_ostream& bgzfStream_ );

    const BamDiffOptions &options_;

    // Input streams
    boost::iostreams::filtering_ostream masterBgzfStream_;
    boost::iostreams::filtering_ostream slaveBgzfStream_;

};


class DiffComputer : public io::bam::BamParserFilter
{
public:
    DiffComputer( unsigned int fileNum );
    ~DiffComputer();

private:
    virtual void parsedRefSeqInfo( const std::vector< BamRefInfoItemType >& bamRefInfo );
    virtual void parsedAlignment( const io::bam::BamAlignment& alignment, const io::bam::VirtualOffset& virtualOffset, const io::bam::VirtualOffset& virtualEndOffset );
    virtual void finishedParsing();

    void ProcessMatch( const io::bam::StorableBamAlignment& masterAlignment, const io::bam::StorableBamAlignment& slaveAlignment );
    void ProcessUnmatchedMaster( const io::bam::StorableBamAlignment& alignment );
    void ProcessUnmatchedSlave( const io::bam::StorableBamAlignment& alignment );
    void UnprocessUnmatchedSlave( const io::bam::StorableBamAlignment& alignment );

    void ProcessCorrect( const io::bam::StorableBamAlignment& masterAlignment, io::bam::StorableBamAlignment& slaveAlignment );
    void ProcessFP( io::bam::StorableBamAlignment& alignment );
    void ProcessFN( const io::bam::StorableBamAlignment& alignment );
    void FinaliseFP();

    void Finalise();

    void ParseReadName( io::bam::StorableBamAlignment& alignment, unsigned int& lane, unsigned int& tile, unsigned long& cluster);
    void ParseReadName( io::bam::StorableBamAlignment& alignment, unsigned int& fullTileNum, unsigned long& cluster);
    void ParseReadNameAndAddSimulatedPosInfo( io::bam::StorableBamAlignment& alignment );


    static std::vector< io::bam::StorableBamAlignment > masterAlignments_;
    static std::vector< io::bam::StorableBamAlignment > slaveAlignments_;
    static std::vector< unsigned int > slave2MasterRefId_;
    static model::FragmentPosResolver fragmentPosResolver_;
    unsigned int fileNum_;
};


} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_BAM_ANALYSER_HH
