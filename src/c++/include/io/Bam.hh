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

#ifndef EAGLE_IO_BAM_HH
#define EAGLE_IO_BAM_HH
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
namespace io
{

// BAM compressor, which we copy-pasted from iSSAC source code and adapted for EAGLE
struct iTag
{
    iTag(const char tag[2], int value):
        val_type_('i'), value_(value)
    {
        tag_[0] = tag[0];
        tag_[1] = tag[1];
    }
    char tag_[2];
    char val_type_;
    int  value_;
}__attribute__ ((packed));

BOOST_STATIC_ASSERT(7 == sizeof(iTag));

struct zTag
{
    zTag(const char tag[2], const char *value):
        val_type_('Z'), value_(value)
    {
        tag_[0] = tag[0];
        tag_[1] = tag[1];
    }
    char tag_[2];
    char val_type_;
    const char *value_;

    size_t size() const {return sizeof(tag_) + sizeof(val_type_) + strlen(value_) + 1;}
}__attribute__ ((packed));

void serialize(std::ostream &os, const char* bytes, size_t size);

inline void serialize(std::ostream &os, const char* pStr) {
    serialize(os, pStr, strlen(pStr) + 1);
}

inline void serialize(std::ostream &os, const std::string &str) {
    serialize(os, str.c_str(), str.length() + 1);
}

//todo: provide proper implementation with byte flipping
inline void serialize(std::ostream &os, const int &i) {
    serialize(os, reinterpret_cast<const char*>(&i), sizeof(i));
}

inline void serialize(std::ostream &os, const char &c) {
    serialize(os, &c, sizeof(c));
}

//todo: provide proper implementation with byte flipping
inline void serialize(std::ostream &os, const unsigned &ui) {
    serialize(os, reinterpret_cast<const char*>(&ui), sizeof(ui));
}

inline void serialize(std::ostream &os, const iTag &tag) {
    serialize(os, tag.tag_, sizeof(tag.tag_));
    serialize(os, tag.val_type_);
    serialize(os, tag.value_);
}

inline void serialize(std::ostream &os, const zTag &tag) {
    serialize(os, tag.tag_, sizeof(tag.tag_));
    serialize(os, tag.val_type_);
    serialize(os, tag.value_);
}


template <typename T>
void serialize(std::ostream &os, const std::vector<T> &vector) {
    serialize(os, reinterpret_cast<const char*>(&vector.front()), vector.size() * sizeof(T));
}

template <typename IteratorT>
void serialize(std::ostream &os, const std::pair<IteratorT, IteratorT> &pairBeginEnd) {
    serialize(os, reinterpret_cast<const char*>(&*pairBeginEnd.first),
              std::distance(pairBeginEnd.first, pairBeginEnd.second) * sizeof(*pairBeginEnd.first));
}

static const unsigned MAX_LANES_PER_FLOWCELL = 8;
static const unsigned MAX_TILES_PER_LANE = 2048;

template <typename THeader>
void serializeHeader(
    std::ostream &os,
    const std::vector<std::string>& argv,
    const THeader &header)
{
    struct Header
    {
        char magic[4];
        int l_text;
    } __attribute__ ((packed));

    const std::string commandLine(boost::join(argv, " "));

    std::string headerText(
        "@HD\t"
        "VN:1.0\t"
        "SO:coordinate\n"
        "@PG\t" 
        "ID:EAGLE\t"
        "PN:EAGLE\t"
        "VN:"
        EAGLE_VERSION
        "\n");
// This was between PN and VN:            "CL:" + commandLine + "\t"
    BOOST_FOREACH(const typename THeader::ReadGroupType &readGroup, header.getReadGroups())
    {
        std::clog << "Writing read group: " << readGroup.getId() << std::endl;
        headerText += readGroup.getValue() + "\n";
    }

    Header bamHeader ={ {'B','A','M',1}, static_cast<int>(headerText.size() + 1) };

    if (!os.write(reinterpret_cast<char*>(&bamHeader), sizeof(bamHeader))){
        BOOST_THROW_EXCEPTION(eagle::common::IoException(errno, "Failed to write BAM header into bam stream"));
    }
    serialize(os, headerText);



//    struct RefSeq
//    {
//        int l_name;
//        char name[5];
//        int l_ref;
//    } __attribute__ ((packed));
//    RefSeq refSeqPhix = {5, "phix", 1234};
//    RefSeq refSeqEcoli = {5, "Ecol", 4321};
//    os.write(reinterpret_cast<char*>(&refSeqPhix), sizeof(refSeqPhix));
//    os.write(reinterpret_cast<char*>(&refSeqEcoli), sizeof(refSeqEcoli));

    serialize(os, header.getRefSequenceCount()); //n_ref
    BOOST_FOREACH(const typename THeader::RefSeqType &refSeq, header.getRefSequences())
    {
        std::clog << "Writing ref sequence name: " << refSeq.name() << " total bases:" << refSeq.length() << " to bam\n";
        const std::string &name(refSeq.name());
        int l_name(name.length() + 1);
        int l_ref(refSeq.length());
        serialize(os, l_name);
        serialize(os, name);
        serialize(os, l_ref);
    }
/*
*/
}

/*!
  @abstract    Calculate the minimum bin that contains a region [beg,end).
  @param  beg  start of the region, 0-based
  @param  end  end of the region, 0-based
  @return      bin
 */
inline int bam_reg2bin(unsigned beg, unsigned end)
{
    --end;
    if (beg>>14 == end>>14) return 4681 + (beg>>14);
    if (beg>>17 == end>>17) return  585 + (beg>>17);
    if (beg>>20 == end>>20) return   73 + (beg>>20);
    if (beg>>23 == end>>23) return    9 + (beg>>23);
    if (beg>>26 == end>>26) return    1 + (beg>>26);
    return 0;
}

//typedef std::vector<flowcell::TileMetadata> TileMetadataList;

template <typename T>
void serializeAlignment(std::ostream &os, T&alignment)
{
//    std::cerr << "writing: aginment\n";

    const int refID(alignment.refId());
    const int pos(alignment.pos());

    const char *readName = alignment.readName();
    const size_t readNameLength = strlen(readName);
    ASSERT_MSG(0xFF > readNameLength, "Read name length must fit in 8 bit value");

    const unsigned bin_mq_nl(unsigned(bam_reg2bin(pos, pos + alignment.seqLen())) << 16 |
                             unsigned(alignment.mapq()) << 8 |
                             unsigned(readNameLength + 1));

    typedef typename T::CigarBeginEnd CigarBeginEnd;
    const CigarBeginEnd cigarBeginEnd = alignment.cigar();
    const size_t cigarLength = std::distance(cigarBeginEnd.first, cigarBeginEnd.second);
    ASSERT_MSG(0xFFFF >= cigarLength, "Cigar length must fit in 16 bit value");

    unsigned flag_nc(unsigned (alignment.flag()) << 16 | (unsigned (cigarLength) & 0xFFFF));
    const int l_seq(alignment.seqLen());
    const int next_RefID(alignment.nextRefId());
    const int next_pos(alignment.nextPos());
    const int tlen(alignment.tlen());

    const std::vector<unsigned char> &seq = alignment.seq();
    const std::vector<unsigned char> &qual = alignment.qual();

/*
    const iTag fragmentSM = alignment.getFragmentSM();
    const iTag fragmentAS = alignment.getFragmentAS();
    const zTag fragmentRG = alignment.getFragmentRG();
    const iTag fragmentNM = alignment.getFragmentNM();
    const iTag fragmentBC = alignment.getFragmentBC();
*/

    const int block_size(  sizeof(refID)
                         + sizeof(pos)
                         + sizeof(bin_mq_nl)
                         + sizeof(flag_nc)
                         + sizeof(l_seq)
                         + sizeof(next_RefID)
                         + sizeof(next_pos)
                         + sizeof(tlen)
                         + readNameLength + 1
                         + cigarLength * sizeof(unsigned)
                         + seq.size()
                         + qual.size()
/*
                         + sizeof(fragmentSM)
                         + sizeof(fragmentAS)
                         + sizeof(fragmentNM)
                         + sizeof(fragmentBC)
                         + fragmentRG.size()
*/
                         );

//    std::cerr << "writing: aginment block size " << block_size << "\n";
    serialize(os, block_size);
    serialize(os, refID);
    serialize(os, pos);

    serialize(os, bin_mq_nl);
    serialize(os, flag_nc);

    serialize(os, l_seq);
    serialize(os, next_RefID);
    serialize(os, next_pos);
    serialize(os, tlen);

    serialize(os, readName);
    serialize(os, cigarBeginEnd);//4
    serialize(os, seq);  //1
    serialize(os, qual); //2

/*
    serialize(os, fragmentSM);
    serialize(os, fragmentAS);
    serialize(os, fragmentRG);
    serialize(os, fragmentNM);
    serialize(os, fragmentBC);
*/
//    std::cerr << "writing: aginment done\n";

}

void serializeBgzfFooter(std::ostream &os);


} // namespace io
} // namespace eagle

#endif // EAGLE_IO_BAM_HH
