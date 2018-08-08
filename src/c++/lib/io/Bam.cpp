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

#include "io/Bam.hh"


using namespace std;

#define INDEX_FOREACH(index,a,b)                 \
    for(int index = -1; index == -1;)            \
        BOOST_FOREACH(a,b) if(++index,true)


namespace eagle
{
namespace io
{


void serialize(std::ostream &os, const char* bytes, size_t size)
{
//    std::cerr << "writing: " << size << " bytes\n";
    if (!os.write(bytes, size))
    {
        BOOST_THROW_EXCEPTION(
            eagle::common::IoException(errno, (boost::format("Failed to write %d bytes into bam stream") % size).str()));
    }
}

void serializeBgzfFooter(std::ostream &os)
{
    // For some strange reason, samtools wants an empty block at the very end of a compressed bam.
    // They call it 'magic'. Note, last \0 is removed (comparing to the original bgzf.c) because
    // the C++ compiler (rightfully) complains.
    const static char magic[28] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0";
    serialize(os, magic, sizeof(magic));
}


} // namespace genome
} // namespace eagle
