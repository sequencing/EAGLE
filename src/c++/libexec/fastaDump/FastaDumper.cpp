/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Top level component to dump a reference.
 **
 ** \author Mauricio Varea
 **/

#include <iostream>
#include <utility>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "model/Struct.hh"
#include "FastaDumper.hh"


namespace eagle
{
namespace main
{


bool less(std::string s1, std::string s2) {return s1.length() < s2.length();}

void FastaDumper::run()
{
    bool overlapContigBoundary = false;
    unsigned long offset = 0;
    unsigned long pos = location_.pos();

    std::vector< std::string > contigNames = reference_.allContigNames();
    unsigned int maxContigNameLength = std::max_element( contigNames.begin(), contigNames.end(), less )->length();
    unsigned int maxPosition = boost::lexical_cast<std::string>(location_.pos() + size_).size();
    while (offset < size_)
    {
        bool prevOverlap = overlapContigBoundary;
        char base = (location_.global ? reference_.get( location_.pos(), offset, overlapContigBoundary )
                                      : reference_.get( location_, offset, overlapContigBoundary ) );
        pos %= reference_.estimatedLength();
        if( 0 == (offset % CONTIG_WIDTH) )
        {
            if (offset) {std::cout << std::endl;}
            display( eagle::model::Locus( reference_.currentChromosome(), pos),
                     maxContigNameLength + 2,
                     maxPosition + 2 );
        } else if (prevOverlap != overlapContigBoundary) {
            pos = 1;
            std::cout << std::endl;
            display( eagle::model::Locus( reference_.currentChromosome(), pos),
                     maxContigNameLength + 2,
                     maxPosition + 2 );
            std::cout << std::setfill(' ') << std::setw( (offset % CONTIG_WIDTH) ) << " ";
        }
        ++offset;
        ++pos;
        std::cout << base;
    }
    std::cout << std::endl;
}


MetaLocus FastaDumper::seek(std::string position)
{
    MetaLocus tmp;
    try {
        tmp = MetaLocus( boost::lexical_cast<unsigned long>(position) );
    } catch (boost::bad_lexical_cast &) {
        tmp = MetaLocus( position );
    }
    return tmp;
}


} // namespace main
} // namespace eagle
