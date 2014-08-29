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
#include "BamAnalyserOptions.hh"


namespace eagle
{
namespace main
{


class BamAnalyser
{
public:
    BamAnalyser (const BamAnalyserOptions &options);
    void run();

private:
    const BamAnalyserOptions &options_;

    boost::iostreams::filtering_ostream bgzfStream_;
};


class BamReadDumper : public io::bam::BamParserFilter
{
    virtual void parsedAlignment( const io::bam::BamAlignment& alignment, const io::bam::VirtualOffset& virtualOffset, const io::bam::VirtualOffset& virtualEndOffset )
    {
        std::cout << "pos=" << alignment.pos << std::endl;
    }
};

class StatsPerGenomeWindow
{
public:
    StatsPerGenomeWindow( unsigned int windowSize) : windowSize_( windowSize ) {}
    void inc( unsigned long globalPos )
        {
            unsigned int index = globalPos / windowSize_;
            if (index >= data_.size())
            {
                data_.resize( index+1 );
            }
            ++(data_[index]);
        }
    void incRegion( unsigned long globalPos, unsigned int length )
        {
            unsigned int index1 = globalPos / windowSize_;
            unsigned int index2 = (globalPos+length-1) / windowSize_;
            if (index2 >= data_.size())
            {
                data_.resize( index2+1 );
            }

            if (index1 == index2)
            {
                (data_[index1]) += length;
            }
            else
            {
                unsigned int length1 = index2 * windowSize_ - globalPos;
                unsigned int length2 = length - length1;
                (data_[index1]) += length1;
                (data_[index2]) += length2;
            }
        }
    unsigned int get( unsigned int index )
        {
            if (index >= data_.size())
            {
                return 0;
            }
            else
            {
                return data_[index];
            }
        }

    std::vector< unsigned int > data_;
    unsigned int windowSize_;
};

class CountingFifo
{
public:
    CountingFifo() : firstCoveredPos_(0), data_(1000), globalPosOfFirstDataElement_(0) {}

    void incRegion( unsigned long globalPos, unsigned int length )
        {
            if (length > data_.size())
            {
                assert(false); // TODO: data_.resize(...) ...
            }
            assert( (globalPos >= firstCoveredPos_) && 
                    (globalPos+length-1 < firstCoveredPos_+data_.size()) );
            for (unsigned int i=0; i<length; ++i)
            {
                ++(data_[(globalPos+i - globalPosOfFirstDataElement_) % data_.size()]);
            }
        }

    unsigned int get( unsigned int globalPos )
        {
            assert( (globalPos >= firstCoveredPos_) && 
                    (globalPos < firstCoveredPos_+data_.size()) );
            return data_[(globalPos - globalPosOfFirstDataElement_) % data_.size()];
        }

    void forgetBefore( unsigned int globalPos )
        {
            while (firstCoveredPos_ < globalPos)
            {
                data_[(firstCoveredPos_ - globalPosOfFirstDataElement_) % data_.size()] = 0;
                ++firstCoveredPos_;
            }
        }

    unsigned long long firstCoveredPos_;

private:
    std::vector< unsigned int > data_;
    unsigned long long globalPosOfFirstDataElement_;
};

class MetricsComputer : public io::bam::BamParserFilter
{
public:
    MetricsComputer( std::vector<bool>& requestedMetrics, std::vector<bool>& requestedTables );

private:
    virtual void parsedAlignment( const io::bam::BamAlignment& alignment, const io::bam::VirtualOffset& virtualOffset, const io::bam::VirtualOffset& virtualEndOffset );
    virtual void finishedParsing();

    std::vector< bool > requestedMetrics_;
    std::vector< bool > requestedTables_;
//    eagle::genome::MultiFastaReference reference_;
    std::vector< std::vector< unsigned int > > mismatches_; // [refBase][seqBase]
    std::vector< std::vector< unsigned int > > homopolymerDeletions_;  // [homopolymerLength][deletionLength]
    std::vector< std::vector< unsigned int > > homopolymerInsertions_; // [homopolymerLength][insertionLength]
    std::vector< unsigned int > homopolymerCount_; // [homopolymerLength]
    std::vector< unsigned int > mismatchCountPerRead_;
    unsigned long mismatchCount_;
    unsigned long insertionCount_;
    unsigned long deletionCount_;
    unsigned long baseCount_;
    StatsPerGenomeWindow mismatchCountPer10kWindow_;
    StatsPerGenomeWindow baseCountPer10kWindow_;
    CountingFifo observedCoverage_;
    CountingFifo readsAddedForTable5_;
    std::deque< std::pair< unsigned long long, unsigned long long > > table5Data_;
    std::vector< unsigned int > insertSizes_;
};


} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_BAM_ANALYSER_HH
