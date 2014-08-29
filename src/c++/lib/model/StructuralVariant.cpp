/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \decription Top level component that deals with variants.
 **
 ** \author Mauricio Varea
 **/

#include <string>
#include <boost/assert.hpp>
#include <boost/io/ios_state.hpp>

#include "model/StructuralVariant.hh"

namespace eagle
{
namespace model
{


std::string StructuralVariant::getTypeName() const
{
    std::string str("");
    if (hasTranslocation())
    {
        str += "Translocation";
    }
    if (hasInsertion() || hasDeletion())
    {
        if (str.length())
            str += " with ";
        if (hasInDel())
            str += "InDel";
        else if (hasInsertion())
            str += "Insertion";
        else
            str += "Deletion";
    }
    if (hasSNP())
    {
        if (str.length())
            str += " and ";
        str += "SNP";
    }
    return str;
}

bool StructuralVariant::operator==( const StructuralVariant& rhs ) const
{
    ComplexRearrangement L = this->getVariant();
    ComplexRearrangement R = rhs.getVariant();
    return (L.adjacency.first.chr() == R.adjacency.first.chr())
        && (L.adjacency.first.decimalPos() == R.adjacency.first.decimalPos())
        && (L.adjacency.second.chr() == R.adjacency.second.chr())
        && (L.adjacency.second.decimalPos() == R.adjacency.second.decimalPos())
        && (L.adjacency.second.dir == R.adjacency.second.dir)
        && (L.sequence.size() == R.sequence.size())
        && std::equal(L.sequence.begin(),L.sequence.end(),R.sequence.begin());
}


std::ostream& operator<<( std::ostream& os, const StructuralVariant& sv )
{
    return os << sv.getVariant() << "\t" << " *" << sv.getTypeName() << "* ";
}


} // namespace model
} // namespace eagle
