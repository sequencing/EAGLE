/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Top level component that deals with variants.
 **
 ** \author Mauricio Varea
 **/

#include <boost/assert.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "genome/Event.hh"


namespace eagle
{
namespace genome
{


#ifdef DISTRIBUTED_GENOME_MUTATOR
int Event::apply2( eagle::model::Contig& contigOut,
                  const EventIterator& lastPosition,
                  const ReferenceBounds& reference,
                  const eagle::model::Direction direction)
{
/*
    std::cout << "... [processing] " << getStructuralVariant() << std::endl;
    EAGLE_DEBUG(8,"from " << lastPosition->getStructuralVariant() );
    EAGLE_DEBUG(8,"in direction " << direction.str() );
*/
    return 0;
}
#endif //ifdef DISTRIBUTED_GENOME_MUTATOR


int Event::apply( eagle::model::Contig& contigOut,
                  const EventIterator& lastPosition,
                  const ReferenceBounds& reference,
                  const eagle::model::Direction direction)
{
    EAGLE_DEBUG(8,"... [processing] " << getStructuralVariant() );
    EAGLE_DEBUG(8,"from " << lastPosition->getStructuralVariant() );
    EAGLE_DEBUG(8,"in direction " << direction.str() );

    // cached contig
    static std::string lastContigId("");
    static ReferenceIterator contig;
    if (lastContigId != src())
    {
        lastContigId = src();
        contig = std::find( reference.first, reference.second, eagle::model::Contig(lastContigId) );
        if (contig == reference.second)
        {
            std::cerr << "ERROR: contig " << lastContigId << " not found!" << std::endl;
            std::cerr << "List of known contigs:" << std::endl;
            for (eagle::genome::ReferenceIterator contig=reference.first; contig != reference.second; ++contig)
            {
                std::cerr << "    " << contig->name() << " (id=" << contig->id() << ")" << std::endl;
            }
        }
    }

    if ( (lastPosition->outgoing().defined() && incoming().defined() && !lastPosition->outgoing().sameAs(incoming()))
    ||  lastPosition->dest() != src() )
    {
        std::stringstream message;
        message << std::endl;
        message << "*** Could not produce a valid DNA segment going from:" << std::endl;
        message << "***       " << lastPosition->getVariant() << std::endl;
        message << "*** To:" << std::endl;
        message << "***       " << getVariant() << std::endl;
        EAGLE_ERROR( message.str() );
    }
    long pos1 = lastPosition->from(direction) + direction.offset();
    long pos2 = to(direction.inv());
    if (0 == pos1) {assert(0 == pos2);}  // only fake event allowed at position 0, until we implement support for telomeres
    std::vector<char> segment = contig->read(pos1,pos2);
    EAGLE_DEBUG(8, "[assign] " << segment.size() << " bases" );
    EAGLE_DEBUG(8, "[assign] " << std::string(segment.begin(), segment.end()).substr(0,100) << (segment.size()>100?"...":"")  ); // Prints max 100 bases
    unsigned long basesCount = contigOut.append( segment, direction.isRev() );
    EAGLE_DEBUG(5,"Copied " << basesCount << " bases");
    if (segment.size() != basesCount)
    {
        EAGLE_WARNING( (boost::format("Only %lu bases copied (%lu expected) while processing:") % basesCount % segment.size()).str() );
        EAGLE_WARNING_CONT( "         " << this->getStructuralVariant() );
        EAGLE_WARNING_CONT( "   from: " );
        EAGLE_WARNING_CONT( "         " << lastPosition->getStructuralVariant() );
    }

    if (!getVariant().sequence.empty())
    {
        // we only need to rev-compl the BiDir events, as everything else has been rev-compl at loading time
        basesCount += contigOut.append( getVariant().sequence, direction.isRev() && getVariant().adjacency.first.dir.isBiDir() );
        EAGLE_DEBUG(10,"+ " << getVariant().sequence.size() << " bases from ALT field (" << std::string( getVariant().sequence.begin(), getVariant().sequence.end() ) << ")");
    }

    // Update metadata DEST field with the current position on the allele as it was before adding inserted bases
    long positionBeforeInsertion = contigOut.size() - getVariant().sequence.size() + (type==model::variant::SNP?1:0);
    metadata_.addInfoValue( "DEST", (boost::format("%s:%d") % contigOut.id() % positionBeforeInsertion).str() );

    return basesCount;
}


std::ostream& operator<<( std::ostream& os, const Event& e )
{
    os << e.getVariant() << " is ";
    if (e.metadata_.id != ".")
    {
        os << "known " << e.metadata_.id;
    } else {
        os << "novel";
    }
    os << " *" << e.getTypeName() << "* in " << e.allele_;
    return os;
}


} // namespace genome
} // namespace eagle
