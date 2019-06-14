/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Component that deals with variant events.
 **
 ** \author Mauricio Varea
 **/

#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "model/Struct.hh"
#include "model/StructuralVariantType.hh"

namespace eagle
{
namespace model
{
namespace variant
{


std::string parseAlternate(const std::string& alt, Breakend& B1, Breakend& B2)
{
    std::string base(alt);
    std::string locus;
    size_t first = alt.find_first_of("[]");
    size_t second = alt.find_first_of("[]",first+1);
    if ( std::string::npos != alt.find_first_of("[]",second+1) )
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF","*** more than two '[' or ']' symbols in the ALT field. Exactly two are needed! ***"));
    }
    if (first != std::string::npos)
    {
        if (second == std::string::npos)
        {
            BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF","*** only one '[' or ']' found in the ALT field. Exactly two are needed! ***"));
        }
        if (alt[first] != alt[second])
        {
            BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF","*** both '[' and ']' are not allowed in the same ALT field ***"));
        }
        if (first)
        {
            base = alt.substr( 0, first );
            locus = alt.substr( first+1, second-first-1 );
            B1.dir = Direction::FWD;
        } else {
            locus = alt.substr( 1, second-1 );
            base = alt.substr( second+1 );
            B1.dir = Direction::REV;
        }
        Locus L(locus);
        if ( 0 == L.chr().length() || (L.chr().find_first_of("<>") != std::string::npos) )
        {
            BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF","*** assembly files are not yet supported ***"));
        }
        B2 = eagle::model::Breakend( locus, (alt[first] == ']' ? Direction::REV : Direction::FWD), base );
    }
    return base;
}

template <>  ComplexRearrangement initialize(std::string chr, unsigned long pos, std::string ref, std::string alt, unsigned int altGtIndex)
{
    if ( 0 == chr.length() || (chr.find_first_of(" \t\n") != std::string::npos) )
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF", (boost::format("*** not a valid CHR field in entry \"%s\t%d\t%s\t%s\" ***") % chr % pos % ref % alt).str()));
    }
    Breakend bnd1(chr,pos);
    Breakend bnd2(chr,pos);
    std::string seq("");

    if ( 0 == ref.length() || (ref.find_first_of(" \t\n") != std::string::npos) )
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF", (boost::format("*** not a valid REF field in entry \"%s\t%d\t%s\t%s\" ***") % chr % pos % ref % alt).str()));
    }
    if ( 0 == alt.length() || (alt.find_first_of(" \t\n") != std::string::npos) )
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF", (boost::format("*** not a valid ALT field in entry \"%s\t%d\t%s\t%s\" ***") % chr % pos % ref % alt).str()));
    }

    if ("." == ref && "." == alt)
    {
        return ComplexRearrangement(bnd1,bnd2,seq,altGtIndex);
    }

    std::string targetBase = parseAlternate( alt, bnd1, bnd2 );
    if ( (ref.find_first_not_of("abcdghkmnrstuvwy.ABCDGHKMNRSTUVWY") != std::string::npos) )
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF",
                             (boost::format("*** '%s' contains 1 or more invalid base(s) for the REF field in entry \"%s\t%d\t%s\t%s\" ***") % ref % chr % pos % ref % alt).str()));
    }
    if ( (targetBase.find_first_not_of("abcdghkmnrstuvwy.ABCDGHKMNRSTUVWY") != std::string::npos) )
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF",
                             (boost::format("*** '%s' contains 1 or more invalid base(s) for the ALT field in entry \"%s\t%d\t%s\t%s\" ***") % targetBase % chr % pos % ref % alt).str()));
    }
    bool jump = (bnd1 != bnd2);
    unsigned long n1 = ref.length() - 1;
    unsigned long n2 = targetBase.length() - 1;
    // Insertion, Deletion, or something else?
    if( !n1 && (boost::starts_with(targetBase,ref) || boost::ends_with(targetBase,ref)) )
    {   // unitary REF
        bnd1.base = ref;
        bnd2.base = ref;
        seq = bnd1.dir.isRev() ? targetBase.substr(0,n2) : targetBase.substr(1);
        Locus delta(chr,0,true);
        bnd1 += delta;
        if (jump)
        {
            bnd2 -= delta;
        } else {
            bnd2 += delta;
        }
    }
    else if( !n2 && (boost::starts_with(ref,targetBase) || boost::ends_with(ref,targetBase)) )
    {   // unitary ALT
        bnd1.base = targetBase;
        bnd2.base = targetBase;
        Locus delta1( chr, bnd1.dir.offset(true), false );
        Locus delta2( chr, ref.length() - 1, false );
        bnd1 += delta1;
        bnd2 += delta2;
    }
    else if( !n1 && !n2 && targetBase[0] != ref[0] )
    {   // SNP
        bnd1.base = ref;
        bnd2.base = "";
        seq = targetBase;
    }
    else // maybe InDel?
    {
        if (!n1 || !n2)
        {
            EAGLE_WARNING( (boost::format("Unitary REF/ALT not correctly detected... {%s:%d, ref='%s', alt='%s'}") % chr % pos % ref % alt).str() );
        }
        if (targetBase[0] == ref[0])
        {
            bnd1.base = ref[0];
            bnd2.base = targetBase[0];
            if (n1 || n2)
                seq = (n1 && !n2) ? ref.substr(1) : targetBase.substr(1);
        } else if(targetBase[n2] == ref[n1]) {
            bnd1.base = ref[n1];
            bnd2.base = targetBase[n2];
            if (n1 || n2)
                seq = (n1 && !n2) ? ref.substr(0,n1) : targetBase.substr(0,n2);
        } else {
            BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF",
                                 (boost::format("*** Variant structure not understood: {%s:%d, ref='%s', alt='%s'} ***") % chr % pos % ref % alt).str() ));
        }
        if (seq.empty())
        {
            EAGLE_WARNING( (boost::format("Empty sequence on a variant that looked like an InDel... {%s:%d, ref='%s', alt='%s'}") % chr % pos % ref % alt).str() );
        }
        Locus delta1( chr, 0, true );
        Locus delta2( chr, ref.length() - 1, true );
        bnd1 += delta1;
        if (jump)
        {
            bnd2 -= delta1;
        } else {
            bnd2 += delta2;
        }
    }
    return ComplexRearrangement( bnd1, bnd2, seq, altGtIndex );
}

template <>  variant::Type initialize(std::string chr, unsigned long pos, std::string ref, std::string alt, unsigned int altGtIndex)
{
    variant::Type svt=variant::Undefined;
    Breakend bnd1(chr,pos);
    Breakend bnd2(chr,pos);
    std::string targetBase = parseAlternate( alt, bnd1, bnd2 );

    if ( (bnd1.chr() != bnd2.chr())
     || (bnd1.pos() != bnd2.pos()) )
    {   // Translocation?
        svt |= variant::Translocation;
    }
    if ( bnd1.chr() == bnd2.chr()
         && (
             bnd1.pos() > bnd2.pos() + 1
          || bnd2.pos() > bnd1.pos() + 1
            )
       )
    {   // Large deletion
        svt |= variant::DEL;
    }
    if (("." == ref && "." == alt) || (0 == ref.length() && 0 == alt.length()))
    {
        return svt;
    }

    unsigned long n1 = ref.length() - 1;
    unsigned long n2 = targetBase.length() - 1;
    // Insertion, Deletion, or SNP?
    if (n1 < n2 && boost::starts_with(targetBase,ref))
    {
        svt |= variant::INS;
    }
    if (n2)
    {
        if( (!n1) || (n1 && targetBase[0] == ref[0] /*&& targetBase[1] != ref[1]*/)
                  || (n1 && targetBase[n2] == ref[n1] /*&& targetBase[n2-1] != ref[n1-1]*/) )
        {
            svt |= variant::INS;
        }
    }
    if (n1 > n2 && boost::starts_with(ref,targetBase))
    {
        svt |= variant::DEL;
    }
    if (n1)
    {
        if( (!n2) || (n2 && targetBase[0] == ref[0] /*&& targetBase[1] != ref[1]*/)
                  || (n2 && targetBase[n2] == ref[n1] /*&& targetBase[n2-1] != ref[n1-1]*/) )
        {
            svt |= variant::DEL;
        }
    }
    if (!n1 && !n2 && targetBase != ref)
    {
        svt |= variant::SNP;
    }

    if (ref=="" && alt=="")
    {   // Special variant marker for start and end of chromosomes
        return svt;
    }

    if (svt == variant::Undefined)
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF",
                             (boost::format("*** Cannot guess variant type for variant: {%s:%d, ref='%s', alt='%s'} ***") % chr % pos % ref % alt).str() ));
    }
    return svt;
}


} // namespace variant
} // namespace model
} // namespace eagle
