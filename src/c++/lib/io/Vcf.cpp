/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Component to read/write VCF files.
 **
 ** \author Mauricio Varea
 **/

#include <algorithm>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "io/Vcf.hh"
#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "model/SplitString.hh"

using namespace boost::lambda;


namespace eagle
{
namespace io
{


VcfMetadata::VcfMetadata(std::string idn, std::string qty, std::string flt, std::string info, std::string format, std::string data)
    : id(idn)
    , qual(qty)
    , filter(flt)
    , infoFieldValue_(info)
    , formatFieldValue_(format)
    , formatDataFieldValue_(data)
    , infoFieldParsed_(false)
    , formatFieldParsed_(false)
{
    if ( 0 == id.length() || (id.find_first_of(" \t\n") != std::string::npos) )
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF",
                                                             (boost::format("*** \"%s\" is not a valid ID field value in entry \"%s\t%s\t%s\t%s\t%s\t%s\t\" ***") % id % idn % qty % flt % info % format % data).str() ));
    }
    if ( 0 == qual.length() || (qual.find_first_not_of("0123456789.") != std::string::npos) )
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF",
                                                             (boost::format("*** \"%s\" is not a valid QUAL field value in entry \"%s\t%s\t%s\t%s\t%s\t%s\t\" ***") % qual % idn % qty % flt % info % format % data).str() ));
    }
    if ( 0 == filter.length() || (filter.find_first_of(" \t\n") != std::string::npos) )
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF",
                                                             (boost::format("*** \"%s\" is not a valid FILTER field value in entry \"%s\t%s\t%s\t%s\t%s\t%s\t\" ***") % filter % idn % qty % flt % info % format % data).str() ));
    }
}

std::string VcfMetadata::strFormat() const
{
    if (formatFieldParsed_)
    {
        std::vector<std::string> keys;
        std::transform( format_.begin(), format_.end(), back_inserter(keys), bind(&InfoType::value_type::first,_1) );
        return boost::algorithm::join( keys, ":" );
    }
    else
    {
        return formatFieldValue_;
    }
}

void VcfMetadata::lazilyParseInfoField() const
{
    if (infoFieldParsed_) { return; }

    typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
    boost::char_separator<char> sep1(";"), sep2("=");

    Tokenizer tok1(infoFieldValue_, sep1);

    for (Tokenizer::iterator it1 = tok1.begin(); it1 != tok1.end(); ++it1)
    {
        Tokenizer tok2(*it1, sep2);
        std::string key;
        std::vector<std::string> values;
        for (Tokenizer::iterator it2 = tok2.begin(); it2 != tok2.end(); ++it2)
        {
            if (key.empty())
            {
                key = *it2;
            } else {
                boost::split( values, *it2, boost::is_any_of(","));
            }
        }
        info_[key] = values;
    }

    infoFieldValue_.clear(); // Free some RAM
    infoFieldParsed_ = true;
}


void VcfMetadata::lazilyParseFormatField() const
{
    if (formatFieldParsed_) { return; }

    typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
    boost::char_separator<char> sep(":");
    Tokenizer tok1(formatFieldValue_, sep);
    Tokenizer tok2(formatDataFieldValue_, sep);

    if ( std::vector<std::string>(tok1.begin(),tok1.end()).size() != std::vector<std::string>(tok2.begin(),tok2.end()).size() )
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF","*** SAMPLE field does not follow the specified FORMAT ***"));
    }

    std::string key;
    std::vector<std::string> values;
    for (Tokenizer::iterator it1 = tok1.begin(), it2 = tok2.begin();
         it2 != tok2.end(); ++it1,++it2)
    {
        key = *it1;
        boost::split( values, *it2, boost::is_any_of(","));
        format_[key] = values;
    }

    formatFieldValue_.clear(); // Free some RAM
    formatDataFieldValue_.clear();
    formatFieldParsed_ = true;
}


VcfVariant::VcfVariant(std::string chrom, std::string pos, std::string id, std::string ref, std::string alt,
                       std::string qual, std::string filter, std::string info, std::string format, std::string data)
: metadata_(id,qual,filter,info,format,data)
{
    this->clear();
    unsigned int ln = 0;
    while( isdigit(pos[ln]) ) ln++;
    if (pos.length() != ln)
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF",
                             (boost::format("*** Non-numeric chars in second field (i.e. POS) of VCF line. Did not understand '%s' ***") % pos).str() ));
    }
    unsigned long lpos;
    try {
        lpos = boost::lexical_cast<unsigned long>(pos);
    } catch (boost::bad_lexical_cast &) {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("VCF",
                             (boost::format("*** Problem converting second field (i.e. POS) of VCF line. Did not understand '%s' ***") % pos).str() ));
    }

    if (alt.find(',') == std::string::npos)
    {
        this->push_back(eagle::model::StructuralVariant(chrom, lpos, ref, alt, 1));
    }
    else
    {
        std::vector<std::string> tokens;
        boost::split(tokens, alt, boost::is_any_of(","));
        unsigned int altGtIndex = 1;
        for(std::vector<std::string>::const_iterator t = tokens.begin(); t != tokens.end(); ++t)
        {
            this->push_back(eagle::model::StructuralVariant(chrom, lpos, ref, *t, altGtIndex++));
        }
    }
}

VcfVariant::VcfVariant(eagle::model::StructuralVariant sv, VcfMetadata meta)
: metadata_(meta)
{
    this->clear();
    this->push_back(sv);
}

std::ostream& operator<<( std::ostream& os, const VcfVariant& vcf )
{
    for(VcfVariant::const_iterator v=vcf.begin(); v!=vcf.end(); ++v)
    {
        if (vcf.begin() != v)
        {
            os << "\n";
        }
        eagle::model::Breakend ref = v->getVariant().adjacency.first;
        eagle::model::Breakend alt = v->getVariant().adjacency.second;
        unsigned long alt_pos;

        os << ref.chr() << "\t";
        if (v->hasTranslocation() || v->hasInsertion())
        {
            os << ref.pos() << "\t";
            alt_pos = alt.pos() + alt.dir.offset();
        } else {
            os << (ref.pos() - alt.dir.offset()) << "\t";
            alt_pos = alt.posAfter();
        }
        os << vcf.metadata_.id << "\t";
        os << ref.base << "\t";

        std::string bracket;
        eagle::model::ComplexRearrangement cr(v->getVariant());
        if (alt.dir.isRev())
        {
            bracket = std::string(1,(char)VcfVariant::REV);
            if (cr.sequence.size())
            {
                cr.inverse();
                alt.base = std::string( cr.sequence.begin(), cr.sequence.end() ) + alt.base;
            }
        } else {
            bracket = std::string(1,(char)VcfVariant::FWD);
            if (cr.sequence.size())
            {
                alt.base += std::string( cr.sequence.begin(), cr.sequence.end() );
            }
        }

        std::string alt_str("");
        if (ref.dir.isRev())
        {
            alt_str = bracket + (boost::format("%s:%lu") % alt.chr() % alt_pos).str() + bracket + alt.base;
        } else if (ref.dir.isFwd()) {
            alt_str = alt.base + bracket + (boost::format("%s:%lu") % alt.chr() % alt_pos).str() + bracket;
        } else {
            alt_str = alt.base;
        }
        os << alt_str << "\t";
        os << vcf.metadata_.qual << "\t";
        os << vcf.metadata_.filter << "\t";
        os << vcf.metadata_.strInfo();
        if (vcf.metadata_.hasData())
        {
            os << "\t" << vcf.metadata_.strFormat();
            os << "\t" << vcf.metadata_.strData();
        }
    }
    return os;
}


bool VcfReader::getNextVariant(VcfVariant &variant, bool filterSnpsOut, bool filterBeginEndMarkersOut)
{
    std::string line;
    while (!(line=getNextLine()).empty())
    {
        model::SplitString tokens( line, "\t" );
//        std::vector<std::string> tokens;
//        boost::split(tokens, line, boost::is_any_of( "\t" ));

        if (filterSnpsOut && tokens.size() > 4)
        {
            if (tokens[3].size() == 1 && tokens[4].size() == 1) // SNP
            {
                if (filterBeginEndMarkersOut
                    || tokens[3] != "."
                    || tokens[4] != "." )
                {
                    continue;
                }
            }
        }

        // Detect and warn about malformed lines
        switch(tokens.size())
        {
        case 5:
            variant = VcfVariant(tokens[0],tokens[1],tokens[2],tokens[3],tokens[4]);
            break;
        case 6:
            variant = VcfVariant(tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5]);
            break;
        case 7:
            variant = VcfVariant(tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6]);
            break;
        case 8:
            variant = VcfVariant(tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6],tokens[7]);
            break;
        case 9:
            EAGLE_WARNING_IF( pathList_.begin() <= thisPath_ && pathList_.end() > thisPath_,
                              (boost::format("Invalid number of tokens (9) at %s:%lu") % *thisPath_ % lineCount_ ).str() );
            variant = VcfVariant(tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6],tokens[7]);
            break;
        case 10:
            variant = VcfVariant(tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6],tokens[7],tokens[8],tokens[9]);
            break;
        default:
            if (5 > tokens.size())
            {
                // end() (and not begin()) included in range, because it's the next (not the current) one.
                EAGLE_WARNING_IF( pathList_.begin() <= thisPath_ && pathList_.end() > thisPath_,
                                  (boost::format("Only %d tokens in %s:%lu") % tokens.size() % *thisPath_ % lineCount_ ).str() );
                EAGLE_WARNING_CONT( "*** LINE IGNORED ***" );
                continue;
            }
            if (10 < tokens.size())
            {
                static bool warningAlreadyIssued = false;
                if (!warningAlreadyIssued && pathList_.begin() <= thisPath_ && pathList_.end() > thisPath_)
                {
                    EAGLE_WARNING( (boost::format("Tokens after column 10 are not parsed by EAGLE (first occurred at %s:%lu)") % *thisPath_ % lineCount_ ).str() );
                    warningAlreadyIssued = true;
                }
                variant = VcfVariant(tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6],tokens[7],tokens[8],tokens[9]);
            }
        }
        return true;
    }
    return false;
}


void VcfWriter::write(const VcfVariant& vcfVariant)
{
    if ( !(*this << vcfVariant) )
    {
        std::stringstream message;
        message << "Failed to write variant:  " << vcfVariant << std::endl
                << "       into file:  " << *thisPath_ << std::endl;
        BOOST_THROW_EXCEPTION(eagle::common::IoException(errno, message.str()));
    }
    std::ofstream::put('\n');
}

void VcfWriter::writeHeader()
{
    if ( !(*this << vcfHeader_) )
    {
        std::stringstream message;
        message << "Failed to write VCF header" << std::endl
                << "       into file:  " << *thisPath_ << std::endl;
        BOOST_THROW_EXCEPTION(eagle::common::IoException(errno, message.str()));
    }
    std::ofstream::put('\n');
}


} // namespace io
} // namespace eagle
