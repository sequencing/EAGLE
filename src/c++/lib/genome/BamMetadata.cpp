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

#include <boost/iostreams/device/file.hpp>

#include "model/PassFilter.hh"
#include "genome/BamMetadata.hh"
#include "io/Bam.hh"
#include "io/BgzfCompressor.hh"
#include "io/BamIndexer.hh"


using namespace std;

#define INDEX_FOREACH(index,a,b)                 \
    for(int index = -1; index == -1;)            \
        BOOST_FOREACH(a,b) if(++index,true)


namespace eagle
{
namespace genome
{


BamOrMetadataOutput::BamOrMetadataOutput( const boost::filesystem::path outFilename, eagle::io::RunInfo &runInfo, PreferredFastaReader* fastaReference )
    : runInfo_( runInfo )
    , fastaReference_( fastaReference?fastaReference:eagle::genome::SharedFastaReference::get() )
{
    init( outFilename );
}

BamOrMetadataOutput::~BamOrMetadataOutput()
{
    flushReorderedAlignmentsUntilPos( std::numeric_limits<unsigned long>::max() );
    bgzfStream_.pop();
    io::serializeBgzfFooter( bamStream_ );
}

void BamOrMetadataOutput::init( const boost::filesystem::path outFilename )
{
    // BAM creation attempt
    {
        const boost::filesystem::path bamPath( outFilename );
        const int bamGzipLevel = boost::iostreams::gzip::best_speed;
        const std::vector<std::string> argv_;

        boost::iostreams::file_sink bamSink(bamPath.string());
        if (!bamSink.is_open()) {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open output BAM file " + bamPath.string()));
        }
        std::clog << "Creating BAM file: " << bamPath << std::endl;

        bamStream_.push(bamSink);

        bgzfStream_.push(eagle::io::bam::BgzfCompressor(bamGzipLevel));

        { // Add BAM Index
            boost::filesystem::path baiPath( outFilename.string() + ".bai" );
            boost::iostreams::file_sink baiSink(baiPath.string());
            bgzfStream_.push(eagle::io::bam::BamIndexer<boost::iostreams::file_sink>(baiSink));
        }

        bgzfStream_.push(bamStream_);

        io::serializeHeader<EagleBamHeaderAdapter>(
            bgzfStream_
            , argv_
            , EagleBamHeaderAdapter( *fastaReference_ )
            );
    }
}

void BamOrMetadataOutput::add( eagle::genome::ReadClusterWithErrors& readClusterWithErrors )
{
    // example:   HSQ1004_81:3:1206:10700:145121  161     c1.fa   10007   1       70M3D24M6S      =       10151   244     TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACAGGCCTAGGCCTAAGCCTAATACTCCAAAT    CCCFFFFFHHHHHJIJJJJJIJIIIJIHIIIJJJIIIHA2BBB3D2B@(.86@C;=@###########################################    BC:Z:0 SM:i:1   AS:i:0  OC:Z:70M3D24M2D6M
    // paired to: HSQ1004_81:3:1206:10700:145121  81      c1.fa   10151   214     36M2D42M1D22M   =       10007   -244    TAACCCTAACCCGAACCCAAACCCTAACCTAACCCTTCCCTCCCCTAACCCTACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT    ##########################################C@.(=.GF;-7'@??(B?EHGAB0IHFBC:IIIHHGJIHFFAJJIFFCHHFDDDFCBB    BC:Z:0  XD:Z:12T5T17^AA$C2TAA11A24^C$18A2C      SM:i:214        AS:i:0

    static unsigned long fragmentNum = 0;

    bool directionIsForward = false;
    unsigned int readNum = 0;
    unsigned int readNum12 = 1;

    BOOST_FOREACH(const eagle::io::ReadDescription &rd, runInfo_.reads)
    {
        if (!rd.isIndex)
        {
            if (! readClusterWithErrors.eFragment_.structure_.getReadInfo( readNum, directionIsForward ) )
            {
                assert( false );
            }

            const vector< unsigned int > &CIGAR = readClusterWithErrors.getCigar( readNum ); // Calling getCigar first so that the lazy evaluation generates it
            string SEQ = readClusterWithErrors.getNucleotideOrQualitySequenceForRead( readNum, true, !directionIsForward );
            unsigned int usedDnaLength = readClusterWithErrors.getUsedDnaLength( readNum );
            if ( usedDnaLength > readClusterWithErrors.eFragment_.fragment_.fragmentLength_ )
            {
                EAGLE_WARNING( "!!! Template length too short (due to simulated deletions in reads): Some reads are running over the end. I may crash if we run over a chromosome end. Please adjust your template length table." );
            }

            unsigned long startPos = readClusterWithErrors.eFragment_.fragment_.startPos_;
            unsigned long endPos   = readClusterWithErrors.eFragment_.fragment_.startPos_ + readClusterWithErrors.eFragment_.fragment_.fragmentLength_ - usedDnaLength;

            string QNAME = (boost::format("FC:%d") % fragmentNum).str(); //readClusterWithErrors.id_;
            bool isPassingFilter = model::PassFilter::isSequencePassingFilter( SEQ );
            unsigned int FLAG  = 0x3 | (readNum12==1?0x40:0x80) | (directionIsForward?0x20:0x10) | (isPassingFilter?0:0x200);
            unsigned long GlobalPos = directionIsForward?startPos:endPos;
            unsigned int MAPQ = 50;
            unsigned long PNEXT = directionIsForward?endPos:startPos;
            long TLEN = readClusterWithErrors.eFragment_.fragment_.fragmentLength_ * (directionIsForward?1:-1);
            string QUAL  = readClusterWithErrors.getNucleotideOrQualitySequenceForRead( readNum, false, !directionIsForward );
            EagleBamAlignmentAdapter eagleBamAlignmentAdapter( GlobalPos, QNAME, FLAG, MAPQ, 
                                                               directionIsForward?CIGAR:(std::vector<unsigned int>(CIGAR.rbegin(),CIGAR.rend())),
                                                               PNEXT, TLEN, SEQ, QUAL, *fastaReference_ );

            if (directionIsForward)
            {
                flushReorderedAlignmentsUntilPos( GlobalPos );
                io::serializeAlignment( bgzfStream_, eagleBamAlignmentAdapter );
            }
            else
            {
                addToReorderedAlignments( eagleBamAlignmentAdapter );
            }

            ++readNum12;
        }
        ++readNum;
    }
    ++fragmentNum;
}


void BamOrMetadataOutput::addRebased( eagle::genome::ReadClusterWithErrors& readClusterWithErrors, const signed long globalPosShift, const unsigned long firstPosToProcess, const unsigned long lastPosToProcess, const bool dropLastBase, const genome::RefToSampleSegment& cigarModifierHelper )
{
    bool directionIsForward = false;
    unsigned int readNum = 0;
    unsigned int readNum12 = 1;
    vector< unsigned int > softClippedCIGAR;

    BOOST_FOREACH(const eagle::io::ReadDescription &rd, runInfo_.reads)
    {
        if (!rd.isIndex)
        {
            if (!readClusterWithErrors.eFragment_.structure_.getReadInfo( readNum, directionIsForward ))
            {
                assert( false );
            }

            const vector< unsigned int > &CIGAR = readClusterWithErrors.getCigar( readNum, dropLastBase ); // Calling getCigar first so that the lazy evaluation generates it
            string SEQ = readClusterWithErrors.getNucleotideOrQualitySequenceForRead( readNum, true, !directionIsForward, dropLastBase );
            unsigned int usedDnaLength = readClusterWithErrors.getUsedDnaLength( readNum, dropLastBase );
            if ( usedDnaLength > readClusterWithErrors.eFragment_.fragment_.fragmentLength_ )
            {
                EAGLE_WARNING( "Template length too short: Some reads are running over the end. Please adjust your template length table. Bam indexing is expected to fail, and you will need to use 'samtools sort' to re-sort the BAM file." );
            }

            unsigned long startPos1 = readClusterWithErrors.eFragment_.fragment_.startPos_;
            unsigned long startPos2 = readClusterWithErrors.eFragment_.fragment_.startPos_ + readClusterWithErrors.eFragment_.fragment_.fragmentLength_ - usedDnaLength;

            signed long GlobalPos = (directionIsForward?startPos1:startPos2) + globalPosShift;
            if (GlobalPos+usedDnaLength-1 >= (signed long)firstPosToProcess && GlobalPos <= (signed long)lastPosToProcess)
            {
                string QNAME = (boost::format("FC:%d") % readClusterWithErrors.eFragment_.fragment_.fragmentNum_).str(); //readClusterWithErrors.id_;
                bool isPassingFilter = model::PassFilter::isSequencePassingFilter( SEQ );
                unsigned int FLAG  = 0x3 | (readNum12==1?0x40:0x80) | (directionIsForward?0x20:0x10) | (isPassingFilter?0:0x200);
                unsigned int MAPQ = 50;
                unsigned long PNEXT = (directionIsForward?startPos2:startPos1) + globalPosShift;
                long TLEN = readClusterWithErrors.eFragment_.fragment_.fragmentLength_ * (directionIsForward?1:-1);
                string QUAL  = readClusterWithErrors.getNucleotideOrQualitySequenceForRead( readNum, false, !directionIsForward, dropLastBase );

                // Soft-clip the CIGAR string if required
                bool throwThisReadAway = false;
                const vector< unsigned int > &reorderedCIGAR = directionIsForward?CIGAR:(std::vector<unsigned int>(CIGAR.rbegin(),CIGAR.rend()));
                unsigned long globalPosAfterSoftClipping = GlobalPos;
                softClippedCIGAR.clear();
                if (GlobalPos < (signed long)firstPosToProcess)
                {
                    bool success = updateLhsCIGAR( reorderedCIGAR, softClippedCIGAR, readClusterWithErrors, cigarModifierHelper, firstPosToProcess, GlobalPos, FLAG, globalPosAfterSoftClipping );
                    throwThisReadAway = !success;
                }
                if (!throwThisReadAway && GlobalPos+usedDnaLength-1 > (signed long)lastPosToProcess)
                {
                    bool success = updateRhsCIGAR( reorderedCIGAR, softClippedCIGAR, readClusterWithErrors, cigarModifierHelper, lastPosToProcess, GlobalPos, GlobalPos+usedDnaLength-1, FLAG, globalPosAfterSoftClipping );
                    throwThisReadAway = !success;
                }

                if (!throwThisReadAway)
                {
                    const vector< unsigned int > &goodCIGAR = softClippedCIGAR.empty()?reorderedCIGAR:softClippedCIGAR;

                    EagleBamAlignmentAdapter eagleBamAlignmentAdapter( globalPosAfterSoftClipping, QNAME, FLAG, MAPQ, goodCIGAR,
                                                                       PNEXT, TLEN, SEQ, QUAL, *fastaReference_ );

                    if (directionIsForward)
                    {
                        flushReorderedAlignmentsUntilPos( globalPosAfterSoftClipping );
                        io::serializeAlignment( bgzfStream_, eagleBamAlignmentAdapter );
                    }
                    else
                    {
                        if (lastPosToProcess && startPos2+globalPosShift <= lastPosToProcess)
                        {
                            addToReorderedAlignments( eagleBamAlignmentAdapter );
                        }
                    }
                }
            }

            ++readNum12;
        }
        ++readNum;
    }
}


bool BamOrMetadataOutput::updateLhsCIGAR( const vector< unsigned int >& reorderedCIGAR, vector< unsigned int >& softClippedCIGAR, genome::ReadClusterWithErrors& readClusterWithErrors, const genome::RefToSampleSegment& cigarModifierHelper, const unsigned long firstPosToProcess, const signed long GlobalPos, unsigned int& FLAG, unsigned long& globalPosAfterSoftClipping )
{
    bool success = true;
    if (cigarModifierHelper.previousSegmentInGroup_)
    {
        // If this segment is connected to the previous one by an indel, we can throw the segment away as it should have already been processed as part of the previous segment
        // ... except if it was an insertion large enough to include the beginning of this segment
        unsigned int interSegmentIns = cigarModifierHelper.samplePos_ - (cigarModifierHelper.previousSegmentInGroup_->samplePos_ + cigarModifierHelper.previousSegmentInGroup_->segmentLengthWithRefDirection_);
//        unsigned int interSegmentDel = cigarModifierHelper.refPos_ - (cigarModifierHelper.previousSegmentInGroup_->refPos_ + cigarModifierHelper.previousSegmentInGroup_->segmentLengthWithRefDirection_);
        if (firstPosToProcess - GlobalPos > interSegmentIns)
        {
            success = false;
        }
    }
    if (!success == false)
    {
        // Soft-clip lhs bases
        unsigned int clippingLength = firstPosToProcess - GlobalPos;
        softClipCIGAR( reorderedCIGAR, clippingLength, softClippedCIGAR );
        FLAG |= 0x100; // secondary alignment flag bit
        globalPosAfterSoftClipping = firstPosToProcess;
    }
    return success;
}

unsigned int splitAndAppendNCigarEntries( const unsigned int length, const vector< unsigned int > &fromCIGAR, vector< unsigned int > &toCIGAR, vector< unsigned int > &reminderCIGAR, unsigned int& insCount, unsigned int& delCount )
{
    unsigned int cigarOp = 0;
    unsigned int cigarOpCount = 0;
    unsigned int remaining = length;
    insCount = delCount = 0;
    BOOST_FOREACH( const unsigned int& cigarVal, fromCIGAR )
    {
        cigarOpCount = cigarVal >> 4;
        cigarOp = cigarVal & 0xF;

        if (remaining > 0)
        {
            if (cigarOp == 1) // 'I' doesn't affect the soft clipping length
            {
                toCIGAR.push_back( cigarVal );
                insCount += cigarOpCount;
            }
            else
            {
                if (cigarOpCount <= remaining)
                {
                    toCIGAR.push_back( cigarVal );
                    if (cigarOp == 2) // 'D' affects the soft clipping length
                    {
                        delCount += cigarOpCount;
                    }
                    remaining -= cigarOpCount;
                }
                else
                {
                    toCIGAR.push_back( remaining << 4 | cigarOp );
                    if (cigarOp == 2) // 'D'
                    {
                        delCount += remaining;
                    }
                    cigarOpCount -= remaining;
                    reminderCIGAR.push_back( cigarOpCount << 4 | cigarOp );
                    remaining = 0;
                }
            }
        }
        else
        {
            reminderCIGAR.push_back( cigarVal );
        }
    }

    return length - remaining;
}

bool BamOrMetadataOutput::updateRhsCIGAR( const vector< unsigned int >& reorderedCIGAR, vector< unsigned int >& softClippedCIGAR, genome::ReadClusterWithErrors& readClusterWithErrors, const genome::RefToSampleSegment& cigarModifierHelper, const unsigned long lastPosToProcess, const signed long GlobalFirstPos, const signed long GlobalLastPos, unsigned int& FLAG, unsigned long& globalPosAfterSoftClipping )
{
    bool success = true;
    if (!cigarModifierHelper.nextSegmentInGroup_)
    {
        // If this segment is NOT connected to the next one by an indel (i.e. it reaches a translocation breakend), we soft clip
        vector< unsigned int > tmpCIGAR = softClippedCIGAR.empty()?reorderedCIGAR:softClippedCIGAR;
        reverse( tmpCIGAR.begin(), tmpCIGAR.end() );
        softClippedCIGAR.clear();

        // soft-clip rhs bases
        unsigned int clippingLength = GlobalLastPos - lastPosToProcess;
        softClipCIGAR( tmpCIGAR, clippingLength, softClippedCIGAR );

        reverse( softClippedCIGAR.begin(), softClippedCIGAR.end() );
    }
    else
    {
        // Otherwise we adjust the CIGAR string to reflect the indel
        vector< unsigned int > tmpCIGAR1, tmpCIGAR2, tmpCIGAR3, tmpCIGAR4;// = softClippedCIGAR.empty()?reorderedCIGAR:softClippedCIGAR;
//        unsigned int clippingLength = GlobalLastPos - lastPosToProcess;
        unsigned int insCount, delCount;
        unsigned int basesToCopy = lastPosToProcess - GlobalFirstPos + 1;
        tmpCIGAR1.clear();
        tmpCIGAR2.clear();
        splitAndAppendNCigarEntries( basesToCopy, softClippedCIGAR.empty()?reorderedCIGAR:softClippedCIGAR, tmpCIGAR1, tmpCIGAR2, insCount, delCount );
        const genome::RefToSampleSegment* cigarModifierHelper2 = &cigarModifierHelper;
        do
        {
            int interSegmentIns = cigarModifierHelper2->nextSegmentInGroup_->samplePos_ - (cigarModifierHelper2->samplePos_ + cigarModifierHelper2->segmentLengthWithRefDirection_);
            int interSegmentDel = cigarModifierHelper2->nextSegmentInGroup_->refPos_ - (cigarModifierHelper2->refPos_ + cigarModifierHelper2->segmentLengthWithRefDirection_);
            if (interSegmentIns > 0)
            {
                tmpCIGAR3.clear();
                tmpCIGAR4.clear();
                unsigned int n = splitAndAppendNCigarEntries( interSegmentIns, tmpCIGAR2, tmpCIGAR4, tmpCIGAR3, insCount, delCount );
                tmpCIGAR2 = tmpCIGAR3;
                tmpCIGAR1.push_back( ( n + insCount - delCount ) << 4 | 1 ); // nI
            }
            if (interSegmentDel > 0)
            {
                tmpCIGAR1.push_back( interSegmentDel << 4 | 2 ); // nD
            }
            unsigned long nextSegmentLength = labs( cigarModifierHelper2->nextSegmentInGroup_->segmentLengthWithRefDirection_ );
            tmpCIGAR3.clear();
            splitAndAppendNCigarEntries( nextSegmentLength, tmpCIGAR2, tmpCIGAR1, tmpCIGAR3, insCount, delCount );
            tmpCIGAR2 = tmpCIGAR3;
            cigarModifierHelper2 = cigarModifierHelper2->nextSegmentInGroup_;
        }
        while (!tmpCIGAR2.empty());
        softClippedCIGAR = tmpCIGAR1;
    }
    return success;
}

void BamOrMetadataOutput::softClipCIGAR( const vector< unsigned int > &CIGAR, unsigned int clippingLength, vector< unsigned int > &softClippedCIGAR )
{
    unsigned int cigarOp = 0;
    unsigned int cigarOpCount = 0;
    softClippedCIGAR.clear();
    cigarOp = 4; // 'S'
    cigarOpCount = clippingLength;
    softClippedCIGAR.push_back( cigarOpCount << 4 | cigarOp );
    unsigned int clippingLengthRemaining = clippingLength;
    BOOST_FOREACH( const unsigned int& cigarVal, CIGAR )
    {
        cigarOpCount = cigarVal >> 4;
        cigarOp = cigarVal & 0xF;

        if (clippingLengthRemaining > 0)
        {
            if (cigarOp == 1) // 'I'
            {
                clippingLength += cigarOpCount;
                softClippedCIGAR[0] = clippingLength << 4 | (softClippedCIGAR[0] & 0xF);
            }
            else
            {
                if (cigarOpCount <= clippingLengthRemaining)
                {
                    if (cigarOp == 2) // 'D'
                    {
                        clippingLength -= cigarOpCount;
                        softClippedCIGAR[0] = clippingLength << 4 | (softClippedCIGAR[0] & 0xF);
                    }
                    clippingLengthRemaining -= cigarOpCount;
                }
                else
                {
                    if (cigarOp == 2) // 'D'
                    {
                        clippingLength -= clippingLengthRemaining;
                        softClippedCIGAR[0] = clippingLength << 4 | (softClippedCIGAR[0] & 0xF);
                    }
                    cigarOpCount -= clippingLengthRemaining;
                    softClippedCIGAR.push_back( cigarOpCount << 4 | cigarOp );
                    clippingLengthRemaining = 0;
                }
            }
        }
        else
        {
            softClippedCIGAR.push_back( cigarVal );
        }
    }
}

void BamOrMetadataOutput::flushReorderedAlignmentsUntilPos( unsigned long GlobalPos )
{
    while (reorderedAlignmentsWindow_.size() >0 && reorderedAlignmentsWindow_.begin()->GlobalPos_ <= GlobalPos)
    {
        io::serializeAlignment( bgzfStream_, *reorderedAlignmentsWindow_.begin() );
        reorderedAlignmentsWindow_.pop_front();
    }
}

void BamOrMetadataOutput::addToReorderedAlignments( EagleBamAlignmentAdapter &alignment )
{
    list<EagleBamAlignmentAdapter>::reverse_iterator itr;
    for (itr=reorderedAlignmentsWindow_.rbegin(); itr!=reorderedAlignmentsWindow_.rend(); ++itr)
    {
        if (itr->GlobalPos_ <= alignment.GlobalPos_)
        {
            reorderedAlignmentsWindow_.insert( itr.base(), alignment );
            return;
        }
    }
    reorderedAlignmentsWindow_.push_front( alignment );
}


} // namespace genome
} // namespace eagle
