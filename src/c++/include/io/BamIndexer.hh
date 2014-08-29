/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description implements a boost iostreams filter based on BamParserFilter that indexes a BAM input stream:
 ** forwards BAM stream to first output and generates a BAI index stream as second output
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_BAM_BAM_INDEXER_HH
#define EAGLE_BAM_BAM_INDEXER_HH

#include <vector>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "io/BamParserFilter.hh"


namespace eagle
{
namespace io
{
namespace bam
{

namespace bios=boost::iostreams;


template<typename Device>
class BamIndexer : public BamParserFilter
{
public:
    BamIndexer(Device& baiSink);
    BamIndexer(const BamIndexer& that);
    void initStructures();

    virtual void parsedRefSeqInfo( const std::vector< BamRefInfoItemType >& bamRefInfo ) { outputBaiHeader( bamRefInfo.size() ); }
    virtual void parsedEndOfChromosome()                            { outputBaiChromosomeIndex(); }
    virtual void parsedAlignment( const BamAlignment& alignment, const VirtualOffset& virtualOffset, const VirtualOffset& virtualEndOffset );

    template<typename Sink> bool flush( Sink& snk );
    template<typename Sink> void close( Sink& snk );

private:
    void outputBaiHeader( const unsigned int bamRefCount );
    void outputBaiFooter();
    void outputBaiChromosomeIndex();
    void addToBinIndex( const unsigned int bin, const VirtualOffset& virtualOffset, const VirtualOffset& virtualEndOffset );
    void addToLinearIndex( const unsigned int pos, const VirtualOffset& virtualOffset );

    Device baiSink_;

    // Bin index
    unsigned int lastIndexedBin_;
    std::vector< std::vector< VirtualOffsetPair > > binIndex_;

    // Linear index
    std::vector< VirtualOffset > linearIndex_;

    // Stats
    unsigned long bamStatsMapped_, bamStatsNmapped_;
};

template<typename Device>
BamIndexer<Device>::BamIndexer(Device& baiSink)
    : BamParserFilter()
    , baiSink_(baiSink)
    , lastIndexedBin_( 0 )
    , bamStatsMapped_( 0 )
    , bamStatsNmapped_( 0 )
{
    binIndex_.resize( BAM_MAX_BIN );
    initStructures();
}

template<typename Device>
BamIndexer<Device>::BamIndexer(const BamIndexer& that)
    : BamParserFilter(that)
    , baiSink_(that.baiSink_)
    , lastIndexedBin_( that.lastIndexedBin_ )
    , binIndex_( that.binIndex_  )
    , linearIndex_( that.linearIndex_  )
    , bamStatsMapped_( that.bamStatsMapped_ )
    , bamStatsNmapped_( that.bamStatsNmapped_ )
{
    initStructures();
}

template<typename Device>
void BamIndexer<Device>::initStructures()
{
    ASSERT_MSG( binIndex_.size() == BAM_MAX_BIN, "Unexpected number of bins in Bam index" );
    for (unsigned int i=0; i<BAM_MAX_BIN; ++i)
    {
        binIndex_[i].reserve( MAX_CLUSTER_PER_INDEX_BIN );
    }

    linearIndex_.reserve( BAM_MAX_CONTIG_LENGTH / 16384 );
}

template<typename Device>
void BamIndexer<Device>::outputBaiHeader( const unsigned int bamRefCount )
{
    bios::write(baiSink_, "BAI\1", 4);
    bios::write(baiSink_, reinterpret_cast<const char*>(&bamRefCount), 4);
}

template<typename Device>
void BamIndexer<Device>::outputBaiFooter()
{
    // output number of coor-less reads (special samtools field)
    const unsigned long long zero = 0;
    bios::write(baiSink_, reinterpret_cast<const char*>(&zero), 8);
}

template<typename Device>
void BamIndexer<Device>::outputBaiChromosomeIndex()
{
    struct {
        uint32_t binNum;
        uint32_t nClusters;
        uint64_t offBeg, offEnd;
        uint64_t mapped, nmapped;
    } __attribute__ ((packed))
          specialBin = { BAM_MAX_BIN, 2, 0xFFFFFFFFFFFFFFFFull, 0, bamStatsMapped_, bamStatsNmapped_ };

    LOCAL_CERR_DEV_TRACE( "outputBaiChromosomeIndex" );
    const unsigned int nBin = std::count_if( binIndex_.begin(), 
                                             binIndex_.end(),
                                             boost::bind(&std::vector< VirtualOffsetPair >::empty, _1) == false )
        + 1; // Add samtools' special bin to the count
    bios::write(baiSink_, reinterpret_cast<const char*>(&nBin), 4);

    unsigned int i=0;
    BOOST_FOREACH( const std::vector< VirtualOffsetPair >& binIndexEntry, binIndex_ )
    {
        if ( !binIndexEntry.empty() )
        {
            bios::write(baiSink_, reinterpret_cast<const char*>(&i), 4);
            const unsigned int nChunk = binIndexEntry.size();
            bios::write(baiSink_, reinterpret_cast<const char*>(&nChunk), 4);
            bios::write(baiSink_, reinterpret_cast<const char*>(&binIndexEntry[0]), nChunk*16);

            // Fill in samtools' "specialBin" bamStats
            if (specialBin.offBeg > binIndexEntry[0].first.get())
            {
                specialBin.offBeg = binIndexEntry[0].first.get();
            }
            if (specialBin.offEnd < binIndexEntry[nChunk-1].second.get())
            {
                specialBin.offEnd = binIndexEntry[nChunk-1].second.get();
            }
        }
        ++i;
    }

    // Writing special samtools bin
    bios::write(baiSink_, reinterpret_cast<const char*>(&specialBin), sizeof(specialBin));

    // n_intv
    const unsigned int nIntv = linearIndex_.size();
    bios::write(baiSink_, reinterpret_cast<const char*>(&nIntv), 4);
    bios::write(baiSink_, reinterpret_cast<const char*>(&linearIndex_[0]), nIntv*8);


    // reset variables to make them ready to process the next chromosome
    bamStatsMapped_ = bamStatsNmapped_ = 0;
    ASSERT_MSG( binIndex_.size() == BAM_MAX_BIN, "Unexpected number of bins in Bam index" );
    BOOST_FOREACH( std::vector< VirtualOffsetPair >& binIndexEntry, binIndex_)
    {
        binIndexEntry.clear();
    }
    linearIndex_.clear();
    LOCAL_CERR_DEV_TRACE( "outputBaiChromosomeIndex completed" );
}

template<typename Device>
template<typename Sink>
bool BamIndexer<Device>::flush(Sink& snk)
{
    bool r1 = dynamic_cast< BamParserFilter* >(this)->flush( snk ); // Call parent's flush()
    bool r2 = bios::flush( baiSink_ );
    return r1 && r2;
}

template<typename Device>
template<typename Sink>
void BamIndexer<Device>::close( Sink& snk )
{
    dynamic_cast< BamParserFilter* >(this)->close( snk ); // Call parent's close()
    outputBaiFooter();
    bios::close( baiSink_ );
}

template<typename Device>
void BamIndexer<Device>::parsedAlignment( const BamAlignment& alignment, const VirtualOffset& virtualOffset, const VirtualOffset& virtualEndOffset )
{
    const unsigned int bin = alignment.binMqNl >> 16;
    const unsigned int flag = alignment.flagNc >> 16;
    const unsigned int startPos = alignment.pos;
    const unsigned int endPos = alignment.pos + alignment.lSeq - 1;

    if (flag & BAM_FUNMAP)
    {
        // Update bamStats for samtools' special bin
        ++bamStatsNmapped_;
    } else {
        addToBinIndex( bin, virtualOffset, virtualEndOffset );
        addToLinearIndex( startPos, virtualOffset );
        addToLinearIndex( endPos, virtualOffset );

        // Update bamStats for samtools' special bin
        ++bamStatsMapped_;
    }
}

template<typename Device>
void BamIndexer<Device>::addToBinIndex( const unsigned int bin, const VirtualOffset& virtualOffset, const VirtualOffset& virtualEndOffset )
{
    ASSERT_MSG( binIndex_.size() == BAM_MAX_BIN, "Unexpected number of bins in Bam index" );
    ASSERT_MSG( bin < BAM_MAX_BIN, "Invalid bin number in uncompressed BAM" );

    if (binIndex_[bin].empty() || 
        ( bin != lastIndexedBin_ && binIndex_[bin].rbegin()->second.compressedOffset() != virtualOffset.compressedOffset() ))
    {
        binIndex_[bin].push_back( std::make_pair( virtualOffset, virtualEndOffset ) );
        LOCAL_CERR_DEV_TRACE( (boost::format("new bin %d for virtual offset %s") % bin % virtualOffset).str() );
        lastIndexedBin_ = bin;
    }
    else
    {
        binIndex_[bin].rbegin()->second = virtualEndOffset;
        LOCAL_CERR_DEV_TRACE( (boost::format("extending bin %d for virtual offset %s") % bin % virtualOffset).str() );
    }
}

template<typename Device>
void BamIndexer<Device>::addToLinearIndex( const unsigned int pos, const VirtualOffset& virtualOffset )
{
    ASSERT_MSG( pos < BAM_MAX_CONTIG_LENGTH, "Alignment position greater than the maximum allowed by BAM index" );
    const unsigned int linearBin = pos>>14;
    if ( linearIndex_.size() <= linearBin )
    {
        const VirtualOffset& lastValue = linearIndex_.empty()?VirtualOffset():linearIndex_.back();
        while ( linearIndex_.size() <= linearBin )
        {
            linearIndex_.push_back( lastValue );
        }
        linearIndex_[linearBin] = virtualOffset;
    }
}



} // namespace bam
} // namespace io
} // namespace eagle


#endif // EAGLE_BAM_BAM_INDEXER_HH
