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

#include "io/BamParserFilter.hh"


namespace eagle
{
namespace io
{
namespace bam
{

namespace bios=boost::iostreams;


BamParserFilter::BamParserFilter()
    : decompressor_( bios::zlib::default_window_bits, GZIP_INTERNAL_RAM )
    , bgzfParserStage_( BGZF_STAGE_INIT )
    , bgzfParserBytesNeeded_( 0 )
    , bgzfBlockCompressedOffset_( 0 )
    , uncompressedOffsetInBgzfBlock_( 0 )
    , bamParserStage_( BAM_STAGE_INIT )
    , bamParserBytesNeeded_( 0 )
    , lastProcessedRefId_( 0 )
    , exceptionDetected_( 0 )
{
    initStructures();
}

BamParserFilter::BamParserFilter(const BamParserFilter& that)
    : decompressor_( that.decompressor_ )
    , bgzfParserStage_( that.bgzfParserStage_ )
    , bgzfParserBytesNeeded_( that.bgzfParserBytesNeeded_ )
    , bgzfBuf_( that.bgzfBuf_ )
    , bgzfBlockCompressedOffset_( that.bgzfBlockCompressedOffset_ )
    , uncompressedOffsetInBgzfBlock_( that.uncompressedOffsetInBgzfBlock_ )
    , bamParserStage_( that.bamParserStage_ )
    , bamParserBytesNeeded_( that.bamParserBytesNeeded_ )
    , bamParserCurrentVirtualOffset_( that.bamParserCurrentVirtualOffset_ )
    , bamParserCurrentVirtualEndOffset_( that.bamParserCurrentVirtualEndOffset_ )
    , bamParserNextVirtualOffset_( that.bamParserNextVirtualOffset_ )
    , lastProcessedRefId_( that.lastProcessedRefId_ )
    , exceptionDetected_( that.exceptionDetected_ )
{
    ASSERT_MSG( !bgzfParserBytesNeeded_, "Bam parser is not expected to be copied while in progress" );
    initStructures();
}


void BamParserFilter::initStructures()
{
    decompressedBam_.reserve( MAX_UNCOMPRESSED_SIZE );
    bgzfBuf_.reserve( MAX_COMPRESSED_SIZE );
}


BamParserFilter::~BamParserFilter()
{
    if (exceptionDetected_) { return; }
    ASSERT_MSG (bgzfBlockCompressedOffset_ == 0, "Bam Parser filter needs to be closed before being destroyed");
}


void BamParserFilter::parseBgzfStream( const char* inputBlock, const std::streamsize src_size )
{
    // Parse incoming compressed stream
    std::streamsize bytesProcessed( 0 );
    std::streamsize bytesLeft     ( src_size );

    while (bytesLeft >= bgzfParserBytesNeeded_)
    {
        ASSERT_MSG( (bgzfParserBytesNeeded_ > 0) || (bgzfParserStage_ == BGZF_STAGE_INIT),
                          "BGZF parser shouldn't be waiting for 0 bytes" );
        bgzfBuf_.insert( bgzfBuf_.end(), inputBlock + bytesProcessed, inputBlock + bytesProcessed + bgzfParserBytesNeeded_ );
        bytesProcessed += bgzfParserBytesNeeded_;
        bytesLeft -= bgzfParserBytesNeeded_;
        bgzfParserBytesNeeded_ = 0;

        switch (bgzfParserStage_)
        {
        case BGZF_STAGE_INIT:
        {
            ASSERT_MSG( bgzfBuf_.size() == 0, "BGZF parser was not initialised correctly" );
            bgzfParserBytesNeeded_ = 18;
            bgzfParserStage_ = BGZF_STAGE_HEADER;
            startedParsing();
            break;
        }
        case BGZF_STAGE_HEADER:
        {
            ASSERT_MSG( bgzfBuf_.size() == 18, "BGZF header was expected to be 18 bytes long but is different" );

            const unsigned int xLen  = (bgzfBuf_[11] << 8) + bgzfBuf_[10];
            const unsigned int bSize = (bgzfBuf_[17] << 8) + bgzfBuf_[16];
            LOCAL_CERR_DEV_TRACE( (boost::format("XLEN=%d") % xLen).str() );
            LOCAL_CERR_DEV_TRACE( (boost::format("BSIZE=%d") % bSize).str() );

            bgzfParserBytesNeeded_ = bSize - xLen - 19;
            bgzfParserStage_ = BGZF_STAGE_BODY;
            break;
        }
        case BGZF_STAGE_BODY:
        {
            LOCAL_CERR_DEV_TRACE( (boost::format("BGZF body received (pos=%d)") % bgzfBuf_.size()).str() );
            bgzfParserBytesNeeded_ = 8;
            bgzfParserStage_ = BGZF_STAGE_FOOTER;
            break;
        }
        case BGZF_STAGE_FOOTER:
        {
            LOCAL_CERR_DEV_TRACE( (boost::format("footer received (pos=%d)") % bgzfBuf_.size()).str() );
            unsigned int uncompressedSize = *(reinterpret_cast<unsigned int*>(&bgzfBuf_[4]));
            LOCAL_CERR_DEV_TRACE( (boost::format("uncompressedSize=%d") % uncompressedSize).str() );
            uncompressedSize++; // Avoid 'unused variable' warning in case warnings are not activated

            ProcessBgzfBlock();

            bgzfParserBytesNeeded_ = 18;
            bgzfBuf_.clear();
            bgzfParserStage_ = BGZF_STAGE_HEADER;
            break;
        }
        }
    }
    bgzfBuf_.insert( bgzfBuf_.end(), inputBlock + bytesProcessed, inputBlock + bytesProcessed + bytesLeft );
    bgzfParserBytesNeeded_ -= bytesLeft;
    LOCAL_CERR_DEV_TRACE( (boost::format("bgzfParserBytesNeeded_=%d") % bgzfParserBytesNeeded_).str() );
}


void BamParserFilter::ProcessBgzfBlock()
{
    LOCAL_CERR_DEV_TRACE( "Processing bgzf block..." );
    unsigned int lastDecompressedBamSize = decompressedBam_.size();
    bios::back_insert_device< std::vector<char> > decompressorSnk(decompressedBam_);
    bgzfBuf_[3] = '\0'; // reset FLG field of BGZF to discard extra subfields, then skip the extra subfield at bytes 11-17
    decompressor_.write(decompressorSnk, reinterpret_cast<const char*>(&bgzfBuf_[0]) , 10);
    decompressor_.write(decompressorSnk, reinterpret_cast<const char*>(&bgzfBuf_[18]), bgzfBuf_.size()-18);

    unsigned int bgzfCompressedSize = bgzfBuf_.size();
    unsigned int bgzfDecompressedSize = decompressedBam_.size() - lastDecompressedBamSize;
    LOCAL_CERR_DEV_TRACE( (boost::format("BGZF block:   compressed size = %d") % bgzfCompressedSize).str() );
    LOCAL_CERR_DEV_TRACE( (boost::format("BGZF block: decompressed size = %d") % bgzfDecompressedSize).str() );
    bgzfDecompressedSize++; // Avoid 'unused variable' warning in case warnings are not activated

    ParseDecompressedBam();

    bgzfBlockCompressedOffset_ += bgzfCompressedSize;
    LOCAL_CERR_DEV_TRACE( "Finished processing bgzf block." );
}


void BamParserFilter::ParseDecompressedBam()
{
    LOCAL_CERR_DEV_TRACE( "Parsing decompressed bam..." );
    char *bamPtr = &decompressedBam_[0];
    std::streamsize bytesLeft( decompressedBam_.size() );

    while (bytesLeft >= bamParserBytesNeeded_)
    {
        unsigned int bytesToParse = bamParserBytesNeeded_;
        bamParserBytesNeeded_ = 0;

        switch (bamParserStage_)
        {
        case BAM_STAGE_INIT:
        {
            bamParserBytesNeeded_ = 8;
            bamParserStage_ = BAM_STAGE_HEADER;
            break;
        }
        case BAM_STAGE_HEADER:
        {
            LOCAL_CERR_DEV_TRACE( "BAM header received" );
            ASSERT_MSG( bamPtr[0] == 'B', "Corrupted uncompressed BAM header" );
            ASSERT_MSG( bamPtr[1] == 'A', "Corrupted uncompressed BAM header" );
            ASSERT_MSG( bamPtr[2] == 'M', "Corrupted uncompressed BAM header" );
            ASSERT_MSG( bamPtr[3] == '\1', "Corrupted uncompressed BAM header" );

            const unsigned int lText = *reinterpret_cast< unsigned int* >(&bamPtr[4]);
            LOCAL_CERR_DEV_TRACE( (boost::format("l_text=%d") % lText).str() );
            bamParserBytesNeeded_ = lText;
            bamParserStage_ = BAM_STAGE_SAM_HEADER_TEXT;
            break;
        }
        case BAM_STAGE_SAM_HEADER_TEXT:
        {
            LOCAL_CERR_DEV_TRACE( "Sam header text received" );
            bamParserBytesNeeded_ = 4;
            bamParserStage_ = BAM_STAGE_REF_SEQ_NUM;
            break;
        }
        case BAM_STAGE_REF_SEQ_NUM:
        {
            LOCAL_CERR_DEV_TRACE( "Ref seq num received" );
            bamRefCount_ = *reinterpret_cast< unsigned int* >(&bamPtr[0]);
            LOCAL_CERR_DEV_TRACE( (boost::format("bamRefCount_=%d") % bamRefCount_).str() );

            ASSERT_MSG (bamRefCount_ > 0, "Invalid number of sequences in uncompressed BAM (n_ref=0)" );
            bamParserStageLoopLeft_ = bamRefCount_;
            bamRefInfo_.clear();
            bamParserBytesNeeded_ = 4;
            bamParserStage_ = BAM_STAGE_REF_NAME_LENGTH;
            break;
        }
        case BAM_STAGE_REF_NAME_LENGTH:
        {
            LOCAL_CERR_DEV_TRACE( "Ref name length received" );
            const unsigned int lName = *reinterpret_cast< unsigned int* >(&bamPtr[0]);
            LOCAL_CERR_DEV_TRACE( (boost::format("l_name=%d") % lName).str() );

            ASSERT_MSG (lName > 0, "Invalid chromosome name in uncompressed BAM" );
            bamParserBytesNeeded_ = lName + 4;
            bamParserStage_ = BAM_STAGE_REF_SEQ_INFO;
            break;
        }
        case BAM_STAGE_REF_SEQ_INFO:
        {
            LOCAL_CERR_DEV_TRACE( "Ref seq info received" );
            LOCAL_CERR_DEV_TRACE( (boost::format("name=%s") % bamPtr).str() );
            unsigned int lRef = *reinterpret_cast< unsigned int* >(&bamPtr[bytesToParse-4]);
            LOCAL_CERR_DEV_TRACE( (boost::format("l_ref=%d") % lRef).str() );
            bamRefInfo_.push_back( std::make_pair( std::string(bamPtr), lRef ) );
            if (--bamParserStageLoopLeft_ > 0)
            {
                bamParserBytesNeeded_ = 4;
                bamParserStage_ = BAM_STAGE_REF_NAME_LENGTH;
            }
            else
            {
                parsedRefSeqInfo( bamRefInfo_ );

                bamParserBytesNeeded_ = 4;
                bamParserStage_ = BAM_STAGE_ALIGNMENT_BLOCK_SIZE;
            }
            break;
        }
        case BAM_STAGE_ALIGNMENT_BLOCK_SIZE:
        {
            LOCAL_CERR_DEV_TRACE( "Alignment block size received" );
            const unsigned int blockSize = *reinterpret_cast< unsigned int* >(&bamPtr[0]);
            LOCAL_CERR_DEV_TRACE( (boost::format("block_size=%d") % blockSize).str() );

            if (bamParserNextVirtualOffset_.get())
            {
                bamParserCurrentVirtualOffset_ = bamParserNextVirtualOffset_;
            }
            else
            {
                bamParserCurrentVirtualOffset_.set( bgzfBlockCompressedOffset_, uncompressedOffsetInBgzfBlock_);
            }
            bamParserBytesNeeded_ = blockSize;
            bamParserStage_ = BAM_STAGE_ALIGNMENT_DATA;
            break;
        }
        case BAM_STAGE_ALIGNMENT_DATA:
        {
            LOCAL_CERR_DEV_TRACE( "Alignment data received" );
            const BamAlignment &alignment(*reinterpret_cast<BamAlignment*>(&bamPtr[0]));
            const unsigned int flag = alignment.flagNc >> 16;

            if (!(flag & BAM_FUNMAP))
            {
                while (alignment.refId != lastProcessedRefId_)
                {
                    ASSERT_MSG (alignment.refId > lastProcessedRefId_,
                                "Chromosome number in BAM alignment is greater than the number of chromosomes declared in BAM header" );
                    parsedEndOfChromosome();
                    lastProcessedRefId_++;
                }
            }

            bamParserCurrentVirtualEndOffset_.set( bgzfBlockCompressedOffset_, uncompressedOffsetInBgzfBlock_ + bytesToParse);
            if ( bytesLeft == bytesToParse )
            {
                LOCAL_CERR_DEV_TRACE( (boost::format("Virtual offset %s at the end of a compressed block...") % bamParserCurrentVirtualEndOffset_).str() );
                bamParserCurrentVirtualEndOffset_.set( bgzfBlockCompressedOffset_ + bgzfBuf_.size(), 0);
                LOCAL_CERR_DEV_TRACE( (boost::format(" ... changed to %s") % bamParserCurrentVirtualEndOffset_).str() );
            }
            parsedAlignment( alignment, bamParserCurrentVirtualOffset_, bamParserCurrentVirtualEndOffset_ );

            bamParserBytesNeeded_ = 4;
            bamParserStage_ = BAM_STAGE_ALIGNMENT_BLOCK_SIZE;
            break;
        }
        }

        uncompressedOffsetInBgzfBlock_ += bytesToParse;
        bamPtr    += bytesToParse;
        bytesLeft -= bytesToParse;
        if (bytesLeft > 0)
        {
            bamParserNextVirtualOffset_.set( bgzfBlockCompressedOffset_, uncompressedOffsetInBgzfBlock_);
        }
    }
    if (bytesLeft == 0)
    {
        bamParserNextVirtualOffset_.set( 0 );
    }

    decompressedBam_.erase( decompressedBam_.begin(), decompressedBam_.begin()+(bamPtr-&decompressedBam_[0]) );
    uncompressedOffsetInBgzfBlock_ = -decompressedBam_.size();
    LOCAL_CERR_DEV_TRACE( (boost::format("Finished parsing decompressed bam. decompressedBam.size=%d , bytesNeeded=%d") % decompressedBam_.size() % bamParserBytesNeeded_).str() );
}


} // namespace bam
} // namespace io
} // namespace eagle
