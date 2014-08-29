/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Basic memory structures for genomic analysis.
 **
 ** \author Mauricio Varea
 **/

#include <iostream>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/assign/list_of.hpp>

#include "common/Exceptions.hh"
#include "model/Struct.hh"
#include "model/Nucleotides.hh"

namespace eagle
{
namespace model
{


Locus::Locus(std::string obj)
{
    std::vector<std::string> tokens;
    boost::split(tokens, obj, boost::is_any_of(":"));
    chr_=tokens[0];
    switch (tokens.size())
    {
    case 1: pos_=0;
            break;
    case 2: pos_=2 * boost::lexical_cast<unsigned long>(tokens[1]);
            break;
    default:
        BOOST_THROW_EXCEPTION(eagle::common::EagleException(0, (boost::format("The string '%s' is not a valid initialization for Locus") % obj).str()));
    }
}

bool Locus::operator<( const Locus& rhs ) const
{
    int chrComp = chr_.compare(rhs.chr());
    bool result = ( chrComp < 0 )
               || ( chrComp == 0 && (this->decimalPos() < rhs.decimalPos()) );
    return result;
}

bool Locus::operator>( const Locus& rhs ) const
{
    int chrComp = chr_.compare(rhs.chr());
    bool result = ( chrComp > 0 )
               || ( chrComp == 0 && (this->decimalPos() > rhs.decimalPos()) );
    return result;
}


std::ostream& operator<<( std::ostream& os, const Locus& loc )
{
    return os << (boost::format("%s:%.1lf") % loc.chr() % loc.decimalPos()).str();
}

std::string Breakend::direction() const
{
    switch(dir.value_)
    {
      case Direction::FWD: return "-->";
      case Direction::REV: return "<--";
      case Direction::BIDIR: return "<->";
      default:  return "---";
    }
}

std::ostream& operator<<( std::ostream& os, const Breakend& bnd )
{
    os << bnd.base << "(" << dynamic_cast<const Locus&>(bnd) << ")";
    return os;
}


void ComplexRearrangement::inverse()
{
    std::swap( this->adjacency.first, this->adjacency.second );
    //this->adjacency.first.inv();
    //this->adjacency.second.inv();
    this->setDirection(eagle::model::Direction::REV);
    if (!sequence.empty())
    {
        std::transform(sequence.begin(), sequence.end(), sequence.begin(), IUPAC(true));
        std::reverse(sequence.begin(), sequence.end());
    }
}


std::ostream& operator<<( std::ostream& os, const ComplexRearrangement& obj )
{
    std::string seq(obj.sequence.begin(),obj.sequence.end());
    return os << obj.adjacency.first << obj.adjacency.first.direction()
              << "{" << seq << "}"
              << obj.adjacency.second.direction() << obj.adjacency.second;
}


} // namespace model
} // namespace eagle
