/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Bam alignment structure that can be stored in vectors
 ** The main BamAlignment cannot, for legacy reasons
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_BAM_STORABLE_BAM_ALIGNMENT_HH
#define EAGLE_BAM_STORABLE_BAM_ALIGNMENT_HH

#include <vector>
#include "model/Nucleotides.hh"
#include "BamParserFilter.hh"


namespace eagle
{
namespace io
{
namespace bam
{

// Based on BamAlignment structure
struct StorableBamAlignment
{
    unsigned int refId;
    unsigned int pos;
    unsigned int binMqNl;
    unsigned int flagNc;
    unsigned int lSeq;
    unsigned int nextRefId;
    unsigned int nextPos;
    unsigned int tLen;
    std::vector<char> allTheRest2;

    StorableBamAlignment( const BamAlignment& alignment )
        : refId    ( alignment.refId),
          pos      ( alignment.pos ),
          binMqNl  ( alignment.binMqNl ),
          flagNc   ( alignment.flagNc ),
          lSeq     ( alignment.lSeq ),
          nextRefId( alignment.nextRefId ),
          nextPos  ( alignment.nextPos ),
          tLen     ( alignment.tLen )
    {
        unsigned int sizeOfAllTheRest = alignment.getLReadName() + 4*alignment.getNCigarOp() + (alignment.lSeq+1)/2 + alignment.lSeq;
        allTheRest2.resize( sizeOfAllTheRest );
        std::copy( &alignment.allTheRest, (&alignment.allTheRest)+sizeOfAllTheRest, allTheRest2.begin() );
    }

    inline unsigned int         getBin()       const { return binMqNl >> 16; }
    inline unsigned int         getMapQ()      const { return ( binMqNl >> 8 ) & 0xFF; }
    inline unsigned int         getLReadName() const { return binMqNl & 0xFF; }
    inline unsigned int         getFlag()      const { return flagNc >> 16; }
    inline unsigned int         getNCigarOp()  const { return flagNc & 0xFFFF; }
    const char*                 getReadName()  const { return &allTheRest2[0]; }
    const unsigned int*         getCigar()     const { return reinterpret_cast<const unsigned int*> (&allTheRest2[0] + getLReadName()); }
    const unsigned char*        getSeq()       const { return reinterpret_cast<const unsigned char*>(&allTheRest2[0] + getLReadName() + 4*getNCigarOp()); }
    const char*                 getQual()      const { return &allTheRest2[0] + getLReadName() + 4*getNCigarOp() + (lSeq+1)/2; }

    std::string getReadNameAsString() const
    {
        std::string result;
        const unsigned int length = getLReadName();
        const char* str = getReadName();
        result.resize(length);
        std::copy( str, str+length, result.begin() );
        return result;
    }

    std::string getCigarAsString() const
    {
        std::string result;
        const unsigned int length = getNCigarOp();
        const unsigned int *cig = getCigar();
        for (unsigned int i=0; i<length; ++i)
        {
            unsigned int opLen = cig[i] << 4;
            unsigned int op    = cig[i] & 0xF;
            const char* CIGAR_LETTERS = "MIDNSHP=X";
            assert( op < 9 );
            result += (boost::format("%d%c") % opLen % CIGAR_LETTERS[op]).str();
        }
        return result;
    }

    std::string getSeqAsString() const
    {
        std::string result;
        model::IUPAC baseConverter;
        const unsigned char *seq = getSeq();
        unsigned int seqLength = lSeq;
        for (unsigned int i=0; i<seqLength; ++i)
        {
            char binBase = (i%2==0)?( seq[i/2] >> 4 ):( seq[i/2] & 0xF );
            char base = baseConverter.binToIupac( binBase );
            result += base;
        }
        return result;
    }

    friend std::ostream& operator<<( std::ostream& os, const StorableBamAlignment& obj )
    {
        return os << (boost::format("{ refId=%d, pos=%d, nextRefId=%d, nextPos=%d, bin=%d, mapq=%d, flag=%d=0x%x, readName=%s, cigar=%s, seq=%s }") % obj.refId % obj.pos % obj.nextRefId % obj.nextPos % obj.getBin() % obj.getMapQ() % obj.getFlag() % obj.getFlag() % obj.getReadNameAsString() % obj.getCigarAsString() % obj.getSeqAsString()).str();
    }

    static int SeqCompare( const io::bam::StorableBamAlignment& p1, const io::bam::StorableBamAlignment& p2 )
    {
        unsigned int lSeq1 = p1.lSeq;
        unsigned int lSeq2 = p2.lSeq;

//        assert( lSeq1 == lSeq2 );
        unsigned int lSeq = lSeq1;
        if (lSeq1 != lSeq2)
        {
            lSeq = std::min<unsigned int>( lSeq1, lSeq2 );
            static bool firstTime = true;
            if (firstTime)
            {
                firstTime = false;
                EAGLE_WARNING( "Comparison of 2 sequences of different sizes: " << lSeq1 << " vs " << lSeq2 << ":\n" <<  p1.getSeqAsString() << "\n vs\n" << p2.getSeqAsString() );
            }
        }

        const unsigned char *seq1 = p1.getSeq();
        const unsigned char *seq2 = p2.getSeq();
        unsigned int nbBytes = (lSeq+1)/2;
        for (unsigned int i=0; i<nbBytes; ++i)
        {
            if (seq1[i] < seq2[i]) { return -1; }
            if (seq1[i] > seq2[i]) { return 1; }
        }
/*
        if (lSeq1 < lSeq2)
        {
            return -1;
        }
        else if (lSeq1 > lSeq2)
        {
            return 1;
        }
*/
        return 0;
    }

    static bool SeqCompareLt( const io::bam::StorableBamAlignment& p1, const io::bam::StorableBamAlignment& p2 )
    {
        return (SeqCompare( p1, p2 ) < 0);
    }

    static int PosCompare( const io::bam::StorableBamAlignment& p1, const io::bam::StorableBamAlignment& p2 )
    {
        if (p1.refId == p2.refId)
        {
            if (p1.pos < p2.pos) { return -1; }
            if (p1.pos > p2.pos) { return 1; }
            return 0;
        }
        else
        {
            return (p1.refId < p2.refId)?-1:1;
        }
    }

    static bool PosCompareLt( const io::bam::StorableBamAlignment& p1, const io::bam::StorableBamAlignment& p2 )
    {
        return (PosCompare( p1, p2 ) < 0);
    }
};


} // namespace bam
} // namespace io
} // namespace eagle


#endif // EAGLE_BAM_STORABLE_BAM_ALIGNMENT_HH
