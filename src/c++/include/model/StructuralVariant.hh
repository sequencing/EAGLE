/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Component that deals with variant events.
 **
 ** \author Mauricio Varea, Lilian Janin
 **/

#ifndef EAGLE_MODEL_STRUCTURAL_VARIANT_HH
#define EAGLE_MODEL_STRUCTURAL_VARIANT_HH

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "model/Struct.hh"
#include "model/StructuralVariantType.hh"

namespace eagle
{
namespace model
{

/*
 * \brief Abstract notion of a Structural Variant
 */
class StructuralVariant
{
public:
    StructuralVariant(ComplexRearrangement cr, variant::Type t=variant::Undefined)
    : variant(cr)
    , type(t)
    {}
    StructuralVariant(std::string chr, unsigned long pos, std::string ref = ".", std::string alt = ".", unsigned int altGtIndex = 1)
        : variant(eagle::model::variant::initialize<ComplexRearrangement>(chr,pos,ref,alt,altGtIndex))
        , type(eagle::model::variant::initialize<variant::Type>(chr,pos,ref,alt))
    {}
    bool operator==( const StructuralVariant& rhs ) const;

    ComplexRearrangement & getVariant() {return variant;}
    variant::Type & getType() {return type;}
    ComplexRearrangement const & getVariant() const {return variant;}
    variant::Type const & getType() const {return type;}
    bool isDefined() const {return type.any();}
    bool hasSNP() const {return (type & variant::SNP) == variant::SNP;}
    bool hasInsertion() const {return (type & variant::INS) == variant::INS;}
    bool hasDeletion() const {return (type & variant::DEL) == variant::DEL;}
    bool hasInDel() const {return (type & variant::InDel) == variant::InDel;}
    bool hasTranslocation() const {return (type & variant::Translocation) == variant::Translocation;}
    bool isBeginEndMarker() const {return type == variant::Undefined;} // TODO: Make this test more precise
    std::string getTypeName() const;

protected:
    ComplexRearrangement variant;
    variant::Type type;
};

std::ostream& operator<<( std::ostream& os, const StructuralVariant& sv );


} // namespace model
} // namespace eagle

#endif // EAGLE_MODEL_STRUCTURAL_VARIANT_HH
