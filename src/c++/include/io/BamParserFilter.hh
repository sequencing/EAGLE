/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Implements a boost iostreams filter that parses a BAM input stream,
 ** forwards it to its output and calls some virtual functions for each BAM item
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_BAM_BAM_PARSER_FILTER_HH
#define EAGLE_BAM_BAM_PARSER_FILTER_HH

#include <vector>
#include <boost/format.hpp>
#include <boost/iostreams/filter/gzip.hpp>

//#include "common/Debug.hh"
#define LOCAL_CERR_DEV_TRACE(x) //std::cerr << x << std::endl;
#define ASSERT_MSG(x,y) assert((x) && y)


namespace eagle
{
namespace io
{
namespace bam
{

namespace bios=boost::iostreams;


// See page 9 of SAM-1.3.pdf file from samtools website for a description of these fields
struct BamAlignment
{
    unsigned int refId;
    unsigned int pos;
    unsigned int binMqNl;
    unsigned int flagNc;
    unsigned int lSeq;
    unsigned int nextRefId;
    unsigned int nextPos;
    unsigned int tLen;
    char allTheRest;

    inline unsigned int         getBin()       const { return binMqNl >> 16; }
    inline unsigned int         getMapQ()      const { return ( binMqNl >> 8 ) & 0xFF; }
    inline unsigned int         getLReadName() const { return binMqNl & 0xFF; }
    inline unsigned int         getFlag()      const { return flagNc >> 16; }
    inline unsigned int         getNCigarOp()  const { return flagNc & 0xFFFF; }
    inline const char*          getReadName()  const { return &allTheRest; }
    inline const unsigned int*  getCigar()     const { return reinterpret_cast<const unsigned int*> (&allTheRest + getLReadName()); }
    inline const unsigned char* getSeq()       const { return reinterpret_cast<const unsigned char*>(&allTheRest + getLReadName() + 4*getNCigarOp()); }
    inline const char*          getQual()      const { return &allTheRest + getLReadName() + 4*getNCigarOp() + (lSeq+1)/2; }
} __attribute__ ((packed));


class VirtualOffset
{
    unsigned long long val_;

public:
    VirtualOffset() : val_(0) {}
    void set( unsigned long long cOffset, unsigned int uOffset) { val_ = (cOffset << 16) | uOffset; }
    void set( unsigned long long val)                           { val_ = val; }
    unsigned long long get()                const               { return val_; }
    unsigned long long compressedOffset()   const               { return val_>>16; }
    unsigned int       uncompressedOffset() const               { return val_ & 0xFFFF; }

    friend std::ostream& operator<<( std::ostream& os, const VirtualOffset& virtualOffset )
    {
        return os << (boost::format("{%d, %d}") % (virtualOffset.val_ >> 16) % (virtualOffset.val_ & 0xFFFF)).str();
    }
};
typedef std::pair< VirtualOffset, VirtualOffset > VirtualOffsetPair;


class BamParserFilter
{
    // 512 Mbases is the longest chromosome length allowed in a BAM index
    static const unsigned int BAM_MAX_CONTIG_LENGTH     = 512*1024*1024; 

    // =(8^6-1)/7+1, as defined in samtools
    static const unsigned int BAM_MAX_BIN               = 37450;         

    // Each non-leaf bin contains 8 sub-bins => we expect a maximum of 7 clusters per bin, but we may sometimes
    // get unlucky and a cluster may be split in two if some reads alternate between 2 bins just when they also
    // reach the end of a BGZF block
    static const unsigned int MAX_CLUSTER_PER_INDEX_BIN = 16;            

    // Single uncompressed BGZF chunks cannot contain more than 65535 bytes. Our uncompressed buffer contains
    // 1 uncompressed BGZF chunk plus the remainder of the previous BGZF chunk = 2 chunks in the worst case
    static const unsigned int MAX_UNCOMPRESSED_SIZE     = 65536*2;       

    // Single compressed BGZF chunk plus the remainder of the previous BGZF chunk = 2 chunks in the worst case.
    // Each compressed chunk is believed to be no larger than its uncompressed data, which is limited to 64KB
    static const unsigned int MAX_COMPRESSED_SIZE       = 65536*2;       

    // Buffer size for boost::gzip_decompressor, including some overhead for its internals
    static const unsigned int GZIP_INTERNAL_RAM = MAX_UNCOMPRESSED_SIZE * 2;

    // BAM format constant
    static const unsigned int BAM_FUNMAP = 4; 

public:
    typedef char char_type;
    struct category : bios::multichar_output_filter_tag, bios::closable_tag, bios::flushable_tag {};
    typedef std::pair< std::string, unsigned long > BamRefInfoItemType;

    BamParserFilter();
    BamParserFilter(const BamParserFilter& that);
    void initStructures();
    virtual ~BamParserFilter();

    template<typename Sink> std::streamsize write(Sink &snk, const char* s, std::streamsize n);
    template<typename Sink> bool flush(Sink& snk);
    template<typename Sink> void close(Sink& snk);

protected:
    void parseBgzfStream( const char* s, const std::streamsize src_size );
    void ProcessBgzfBlock();
    void ParseDecompressedBam();

    virtual void startedParsing() {}
    virtual void parsedRefSeqInfo( const std::vector< BamRefInfoItemType >& bamRefInfo ) {}
    virtual void parsedEndOfChromosome() {}
    virtual void parsedAlignment( const BamAlignment& alignment, const VirtualOffset& virtualOffset, const VirtualOffset& virtualEndOffset ) {}
    virtual void finishedParsing() {}

private:
    bios::gzip_decompressor decompressor_;

    // BGZF parser
    enum {BGZF_STAGE_INIT, BGZF_STAGE_HEADER, BGZF_STAGE_BODY, BGZF_STAGE_FOOTER} bgzfParserStage_;
    unsigned int bgzfParserBytesNeeded_;
    std::vector<unsigned char> bgzfBuf_;
    unsigned long long bgzfBlockCompressedOffset_;
    unsigned int uncompressedOffsetInBgzfBlock_;

    // From BGZF to BAM parsers
    std::vector<char> decompressedBam_;

    // BAM parser
    enum {
        BAM_STAGE_INIT
        , BAM_STAGE_HEADER
        , BAM_STAGE_SAM_HEADER_TEXT
        , BAM_STAGE_REF_SEQ_NUM
        , BAM_STAGE_REF_NAME_LENGTH
        , BAM_STAGE_REF_SEQ_INFO
        , BAM_STAGE_ALIGNMENT_BLOCK_SIZE
        , BAM_STAGE_ALIGNMENT_DATA
    } bamParserStage_;
    unsigned int bamParserBytesNeeded_;
    unsigned int bamParserStageLoopLeft_;
    VirtualOffset bamParserCurrentVirtualOffset_;
    VirtualOffset bamParserCurrentVirtualEndOffset_;
    VirtualOffset bamParserNextVirtualOffset_;
    unsigned int bamRefCount_;
    std::vector< BamRefInfoItemType > bamRefInfo_;
    unsigned int lastProcessedRefId_;
    bool exceptionDetected_;
};


template<typename Sink>
std::streamsize BamParserFilter::write(Sink &snk, const char* s, std::streamsize src_size)
{
    if (exceptionDetected_) { return 0; } // Needed because this function will be called by the final stream flush before destruction in case of an Exception

    // Try-catch block at this level to catch exceptions thrown in any of the user callbacks (or in this code)
    try
    {
        ASSERT_MSG( &snk != 0, "Sink missing after Bam Parser filter" );
        LOCAL_CERR_DEV_TRACE( (boost::format("Writing %d bytes to BamParserFilter") % src_size).str() );
        std::streamsize bytesWritten = bios::write(snk, s, src_size);
        ASSERT_MSG( bytesWritten == src_size, "Could not transfer all bytes from BAM source to BAM output" );
        bytesWritten++; // to avoid the stupid 1-page-long warning
        parseBgzfStream( s, src_size );
    }
    catch (...)
    {
        exceptionDetected_ = true;
        throw;
    }

    return src_size;
}

template<typename Sink>
bool BamParserFilter::flush(Sink& snk)
{
    ASSERT_MSG( &snk != 0, "Sink missing after Bam Parser filter" );
    bool ret = bios::flush( snk );
    return ret;
}

template<typename Sink> void BamParserFilter::close( Sink& )
{
    if (exceptionDetected_) { return; }

    if (bgzfBlockCompressedOffset_)
    {
        while (lastProcessedRefId_ != bamRefCount_)
        {
            LOCAL_CERR_DEV_TRACE( (boost::format("lastProcessedRefId_=%d") % lastProcessedRefId_).str() );
            LOCAL_CERR_DEV_TRACE( (boost::format("bamRefCount_=%d") % bamRefCount_).str() );
            ASSERT_MSG (lastProcessedRefId_ < bamRefCount_,
                        "Bam indexer processed more chromosomes than was declared in Bam header" );
            parsedEndOfChromosome();
            lastProcessedRefId_++;
        }
        bgzfBlockCompressedOffset_ = 0;
        finishedParsing();
    }

    // Do I need to propagate the close() to the next Sink?
}



} // namespace bam
} // namespace io
} // namespace eagle


#endif // EAGLE_BAM_BAM_PARSER_FILTER_HH
