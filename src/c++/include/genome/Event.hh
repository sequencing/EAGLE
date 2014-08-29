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

#ifndef EAGLE_GENOME_EVENT_HH
#define EAGLE_GENOME_EVENT_HH

#include <time.h>
#include <string>
#include <utility>
#include <boost/format.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/ptr_container/ptr_deque.hpp>

#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "io/Vcf.hh"
#include "model/Contig.hh"
#include "model/Struct.hh"
#include "model/StructuralVariant.hh"
#include "model/Genotype.hh"
#include "genome/Reference.hh"

namespace eagle
{
namespace genome
{

typedef boost::numeric::interval_lib::rounded_math< unsigned long > rounding_policy;
typedef boost::numeric::interval_lib::checking_no_nan< unsigned long > checking_policy;
typedef boost::numeric::interval_lib::policies< rounding_policy, checking_policy > interval_policies;
typedef boost::numeric::interval< unsigned long, interval_policies > interval;

class Event;
typedef std::vector< Event >::iterator  EventIterator;


class Event : public eagle::model::StructuralVariant
{
  public:
    Event(eagle::model::StructuralVariant sv, eagle::io::VcfMetadata meta, unsigned int ploidy=1)
        : eagle::model::StructuralVariant(sv), metadata_(meta), allele_(ploidy), pairedEvent_(0) {}
    Event(eagle::model::StructuralVariant sv, unsigned int ploidy)
        : eagle::model::StructuralVariant(sv), allele_(ploidy), pairedEvent_(0) {}
    Event(eagle::model::StructuralVariant sv)
        : eagle::model::StructuralVariant(sv), pairedEvent_(0) {}
    Event()
        : eagle::model::StructuralVariant("",0,0), pairedEvent_(0) {}
    bool operator<( const Event & rhs ) const { return (getVariant().adjacency.first < rhs.getVariant().adjacency.first); }
    bool operator==( const Event & rhs ) const { return (getStructuralVariant() == rhs.getStructuralVariant()); }
    bool operator>( const Event & rhs ) const { return (getVariant().adjacency.first > rhs.getVariant().adjacency.first); }
    static bool ltComparisonIncludingAltField( const Event & lhs, const Event & rhs ) {
        return (lhs.getVariant().adjacency.first < rhs.getVariant().adjacency.first) ||
            ( (lhs.getVariant().adjacency.first == rhs.getVariant().adjacency.first) &&
              (lhs.getVariant().sequence < rhs.getVariant().sequence) );
    }
    eagle::model::StructuralVariant getStructuralVariant() const
    {
        return eagle::model::StructuralVariant(this->getVariant(),this->getType());
    }
    eagle::model::Direction incoming() const {return this->getVariant().adjacency.first.dir;}
    eagle::model::Direction outgoing() const {return this->getVariant().adjacency.second.dir;}
    std::string src() const {return this->getVariant().adjacency.first.chr();}
    std::string dest() const {return this->getVariant().adjacency.second.chr();}
#ifdef DISTRIBUTED_GENOME_MUTATOR
    int apply2( eagle::model::Contig& contigOut,
               const EventIterator& lastPosition,
               const ReferenceBounds& reference,
               const eagle::model::Direction direction);
#endif //ifdef DISTRIBUTED_GENOME_MUTATOR
    int apply( eagle::model::Contig& contigOut,
               const EventIterator& lastPosition,
               const ReferenceBounds& reference,
               const eagle::model::Direction direction);
    unsigned long from(const eagle::model::Direction dir = eagle::model::Direction::NONE)
    {
        if (dir.defined())
        {
            return getVariant().adjacency.second.pos(dir);
        }
        return getVariant().adjacency.second.pos();
    }
    unsigned long to(const eagle::model::Direction dir = eagle::model::Direction::NONE)
    {
        if (dir.defined())
        {
            return getVariant().adjacency.first.pos(dir);
        }
        return getVariant().adjacency.first.pos();
    }

    eagle::io::VcfMetadata metadata_;
    eagle::model::Genotype allele_;
    int pairedEvent_;
};

std::ostream& operator<<( std::ostream& os, const Event& v );


} // namespace genome
} // namespace eagle

#endif // EAGLE_GENOME_EVENT_HH
