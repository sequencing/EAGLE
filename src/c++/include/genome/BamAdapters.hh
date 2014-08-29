/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Adapter components for BAM module
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_GENOME_BAM_ADAPTERS_HH
#define EAGLE_GENOME_BAM_ADAPTERS_HH


#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "io/RunInfo.hh"
#include "io/BamParserFilter.hh"
#include "genome/ReadCluster.hh"
#include "genome/SharedFastaReference.hh"
#include "common/Exceptions.hh"
#include "io/StorableBamAlignment.hh"


namespace eagle
{
namespace genome
{

class EagleBamHeaderAdapter
{
public:
    EagleBamHeaderAdapter( PreferredFastaReader& fastaReference )
    {
        const std::vector<std::string> names = fastaReference.allContigNames();
        const std::vector<unsigned long> lengths = fastaReference.allContigLengths();
        for (unsigned int i=0; i<names.size(); ++i) {
            refSequences.push_back( RefSequence( names[i], lengths[i] ) );
        }
    }

    EagleBamHeaderAdapter( const std::vector< io::bam::BamParserFilter::BamRefInfoItemType >& bamRefInfo )
    {
        for (unsigned int i=0; i<bamRefInfo.size(); ++i) {
            refSequences.push_back( RefSequence( bamRefInfo[i].first, bamRefInfo[i].second ) );
        }
    }

    class RefSequence
    {
    public:
        RefSequence( const std::string& name, const int length ) : name_(name), length_(length) {}
        const std::string &name() const {return name_;}
        int length() const {return length_;}
    private:
        std::string name_;
        int length_;
    };

    int getRefSequenceCount() const {
        return refSequences.size();
    }

    typedef RefSequence RefSeqType;
    std::vector<RefSequence> getRefSequences() const {
        return refSequences;
    }

    struct ReadGroupType : public std::map<std::string, std::string>::value_type
    {
        ReadGroupType(const std::map<std::string, std::string>::value_type &that):
            std::map<std::string, std::string>::value_type(that){}
        const std::string &getId() const {return first;}
        const std::string &getValue() const {return second;}
    };

    std::map<std::string, std::string> getReadGroups() const
    {
        std::map<std::string, std::string> ret;
        return ret;
    }

private:
    std::vector<RefSequence> refSequences;
};


class EagleBamAlignmentAdapter
{
public:
    unsigned long GlobalPos_;
    string QNAME_;
    unsigned int FLAG_;
    int refID_;
    unsigned long POS_;
    unsigned int MAPQ_;
    std::vector<unsigned int> CIGAR_;
    int next_refID_;
    unsigned long PNEXT_;
    long TLEN_;
//    string SEQ_;
//    string QUAL_;

    std::vector<unsigned char> SEQ2_;
    std::vector<unsigned char> QUAL2_;
//    std::vector<unsigned int> CIGAR2_;

public:
    EagleBamAlignmentAdapter(
        unsigned long GlobalPos
        , string QNAME
        , unsigned int FLAG
        , unsigned int MAPQ
        , const std::vector<unsigned int>& CIGAR
        , unsigned long GlobalPNEXT
        , long TLEN
        , string SEQ
        , string QUAL
        , PreferredFastaReader& fastaReference
        )
        : GlobalPos_(GlobalPos)
        , QNAME_(QNAME)
        , FLAG_(FLAG)
        , MAPQ_(MAPQ)
        , CIGAR_(CIGAR)
        , TLEN_(TLEN)
    {
        eagle::model::IUPAC converter;

        unsigned int seqSize = SEQ.size();

        // CIGAR
//        CIGAR2_.push_back( seqSize<<4 | 0 ); // 101 M (M=0)

        // TODO: optimise those SEQ2 and QUAL2 away
        // SEQ
        SEQ2_.resize( (seqSize+1) / 2 );
        for (unsigned int i = 0; i<seqSize; ++i)
        {
            if (i%2 == 0)
            {
                SEQ2_[i/2] = converter.bin( SEQ[i] ) << 4;
            }
            else
            {
                SEQ2_[i/2] |= converter.bin( SEQ[i] );
            }
        }

        // QUAL
        QUAL2_.resize( seqSize );
        for (unsigned int i = 0; i<seqSize; ++i)
        {
            QUAL2_[i] = QUAL[i] - 33;
        }

        // {refID_,POS_}
        fastaReference.convertFromGlobalPos( GlobalPos, refID_, POS_);
        --POS_;

        // {next_refID_,PNEXT_}
        fastaReference.convertFromGlobalPos( GlobalPNEXT, next_refID_, PNEXT_);
        --PNEXT_;
    }

    EagleBamAlignmentAdapter( const io::bam::StorableBamAlignment& alignment )
        : QNAME_( alignment.getReadNameAsString() )
        , FLAG_( alignment.getFlag() )
        , refID_( alignment.refId )
        , POS_( alignment.pos )
        , MAPQ_( alignment.getMapQ() )
        , next_refID_( alignment.nextRefId )
        , PNEXT_( alignment.nextPos )
        , TLEN_( alignment.tLen )
    {
        unsigned int seqSize = alignment.lSeq;

        const unsigned char *seq = alignment.getSeq();
        SEQ2_.resize( (seqSize+1) / 2 );
        std::copy( seq, seq+SEQ2_.size(), SEQ2_.begin() );

        const char* qual = alignment.getQual();
        QUAL2_.resize( seqSize );
        std::copy( qual, qual+QUAL2_.size(), QUAL2_.begin() );

        const unsigned int *cig = alignment.getCigar();
        const unsigned int cigLength = alignment.getNCigarOp();
        CIGAR_.resize( cigLength );
        std::copy( cig, cig+cigLength, CIGAR_.begin() );
    }

    static size_t getMaxReadNameLength()
    {
        // TODO: get this calculated and controlled
        return 1024;
    }

    //TODO: get the properly formatted qname
    const char *readName() {
        return QNAME_.c_str();
    }

    typedef std::pair<const unsigned *, const unsigned *> CigarBeginEnd;

    CigarBeginEnd cigar() const {
        return std::make_pair((unsigned*)&CIGAR_[0], (unsigned*)&CIGAR_[CIGAR_.size()]);
    }

    int seqLen() const {
        return QUAL2_.size();
    }

    static unsigned char bamBaseFromBclByte(unsigned char bclByte){
        return bclByte ? 1 << (bclByte & 0x03) : 15;
    }

    static unsigned char bamBasesFromBclShort(unsigned short bclShort){
        unsigned char *pBclShort(reinterpret_cast<unsigned char*>(&bclShort));
        return bamBaseFromBclByte(pBclShort[0]) << 4 | bamBaseFromBclByte(pBclShort[1]);
    }

    static unsigned char bamQualFromBclByte(unsigned char bclByte){
        return bclByte >> 2;
    }

    const std::vector<unsigned char> &seq() {
        return SEQ2_;
    }

    const std::vector<unsigned char> &qual() {
        return QUAL2_;
    }

    int refId() const {
        return refID_;
    }

    int pos() const {
        return POS_;
    }

    unsigned char mapq() const {
        return MAPQ_;
    }

    unsigned flag() const {
        return FLAG_;
    }
    int nextRefId() const {
        return next_refID_;
    }

    int nextPos() const {
        return PNEXT_;
    }
    // todo: do the template length
    int tlen() const {
        return TLEN_;
    }

};


} // namespace genome
} // namespace eagle

#endif // EAGLE_GENOME_BAM_ADAPTERS_HH
