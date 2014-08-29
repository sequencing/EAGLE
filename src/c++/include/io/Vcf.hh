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

#ifndef EAGLE_IO_VCF_HH
#define EAGLE_IO_VCF_HH

#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>


#include "io/Text.hh"
#include "model/StructuralVariant.hh"

namespace eagle
{
namespace io
{

// TODO:     column_names = [ 'Chrom', 'Pos', 'ID', 'Ref', 'Alt', 'Qual', 'Filter', 'Info', 'Format', 'data[0]', 'data[1]' ... ]
// TODO:                    ( only data[0] supported at the moment, i.e. mono-sample )

class VcfMetadata
{
private:
    void lazilyParseInfoField() const;
    void lazilyParseFormatField() const;

public:
    VcfMetadata(std::string idn=".",std::string qty=".",std::string flt="PASS",std::string info="",std::string format="",std::string data="");
    std::string id;
    std::string qual;
    std::string filter;

    std::string strInfo() const { lazilyParseInfoField(); return get<';',true>(info_);}
    std::string strFormat() const;
    std::string strData() const { lazilyParseFormatField(); return get<':',false>(format_);}
    bool hasInfo() const { return !info_.empty() || !infoFieldValue_.empty();}
    bool hasData() const { return !format_.empty() || !formatFieldValue_.empty();}
    void addInfoValue(const std::string &key, const std::string &value)  { lazilyParseInfoField(); info_[key].push_back(value);}
    void addFormatData(const std::string &key, const std::string &value) { lazilyParseFormatField(); format_[key].push_back(value);}
    std::vector< std::string > getInfo(const std::string &key) { lazilyParseInfoField(); return info_[key];}
    std::vector< std::string > getData(const std::string &key) { lazilyParseFormatField(); return format_[key];}

private:
    typedef std::map< std::string, std::vector<std::string> > InfoType;
    template <char C,bool B>
             std::string get(const InfoType& info) const
    {
        std::vector<std::string> defs;
        BOOST_FOREACH(const InfoType::value_type& infoItem, info)
        {
            if (!infoItem.second.empty())
            {
                std::vector<std::string> values(infoItem.second);
                std::stable_sort(values.begin(),values.end());                                // sort values
                values.resize( std::unique(values.begin(),values.end()) - values.begin() );   // remove duplicates
                std::string prefix = infoItem.first + "=";
                defs.push_back( (B ? prefix : "") + boost::algorithm::join( values, "," ) );
            } else {
                if (B) defs.push_back( infoItem.first );
            }
        }
        return boost::algorithm::join( defs, std::string(1,C) );
    }

    // lazy evaluation of fields
    mutable std::string infoFieldValue_;
    mutable std::string formatFieldValue_;
    mutable std::string formatDataFieldValue_;

    mutable InfoType info_;
    mutable InfoType format_;

    mutable bool infoFieldParsed_;
    mutable bool formatFieldParsed_;
};


class VcfVariant : public std::vector<eagle::model::StructuralVariant>
{
public:
    enum Direction {NONE=0,
         FWD='[', REV=']'};
    VcfVariant(std::string chrom, std::string pos, std::string id, std::string ref, std::string alt,
               std::string qual=".", std::string filter="PASS", std::string info="", std::string format="", std::string data="");
    VcfVariant(eagle::model::StructuralVariant sv, VcfMetadata meta);
    VcfVariant() {}

    VcfMetadata metadata_;
};

std::ostream& operator<<( std::ostream& os, const VcfVariant& vcf );


class VcfReader : public DsvReader
{
public:
    VcfReader(const std::vector<boost::filesystem::path> &vcfPathList)
    : DsvReader(vcfPathList) {}
    VcfReader() {}
    bool getNextVariant(VcfVariant &variant, bool filterSnpsOut=false, bool filterBeginEndMarkersOut = false);
};


class VcfWriter : public DsvWriter
{
public:
    VcfWriter(const std::vector<boost::filesystem::path> vcfPathList, const bool overwrite)
    : DsvWriter(vcfPathList,overwrite), vcfHeader_("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILT\tINFO\tFORMAT\tSAMPLE") {}
    VcfWriter(const boost::filesystem::path vcfPath, const bool overwrite)
    : DsvWriter(vcfPath,overwrite),     vcfHeader_("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILT\tINFO\tFORMAT\tSAMPLE") {}
    VcfWriter(const bool overwrite = false)
    : DsvWriter(overwrite),             vcfHeader_("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILT\tINFO\tFORMAT\tSAMPLE") {}
    void write(const VcfVariant& variant);
    void writeHeader();
private:
    std::string vcfHeader_;  // Hardcoded, at the moment (TODO: take this from input file + incorporate changes)
};


} // namespace io
} // namespace eagle

#endif // EAGLE_IO_VCF_HH
