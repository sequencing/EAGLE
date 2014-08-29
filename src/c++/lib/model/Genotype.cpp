/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Syntactic sugar for Genotype info.
 **
 ** \author Mauricio Varea
 **/

#include <string>
#include <algorithm>
#include <boost/assert.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "common/Exceptions.hh"
#include "model/Genotype.hh"

namespace eagle
{
namespace model
{


unsigned int Ploidy::level(const std::string& chrom) const
{
    unsigned int n(organismPloidy_);
    for (std::map<std::string,unsigned int>::const_iterator it = nonPloidyChromosomes_.begin();
         it != nonPloidyChromosomes_.end();
         ++it)
    {
        if (chrom == it->first)
        {
            n = it->second;
            break;
        }
    }
    return n;
}

unsigned int Ploidy::max() const
{
    unsigned int n(organismPloidy_);
    for (std::map<std::string,unsigned int>::const_iterator it = nonPloidyChromosomes_.begin();
         it != nonPloidyChromosomes_.end();
         ++it)
    {
        if (n < it->second)
        {
            n = it->second;
        }
    }
    return n;
}

std::string Ploidy::label(unsigned int p) const
{
    switch(p)
    {
      case 1: return "haploid";
      case 2: return "diploid";
      case 3: return "triploid";
      case 4: return "tetraploid";
      case 5: return "petaploid";
      case 6: return "hexaploid";
     default: return (boost::format("%d-ploid") % p).str();
    }
}

std::ostream& Genotype::dump(std::ostream& os)
{
    for(std::set<int>::iterator it = this->begin(); it != this->end(); ++it)
    {
        os << "\t" << *it << std::endl;
    }
    return os;
}

void Genotype::setPloidy(unsigned int n)
{
    if (n < altSize())
    {
        std::stringstream message;
        message << (boost::format("*** Genotype cannot be set to '%d'. This event already contains %d alleles ***") % n % altSize()).str() << std::endl;
        message << "    The following alleles have already been used:" << std::endl;
        message << dump(message);
        BOOST_THROW_EXCEPTION(common::PreConditionException( message.str() ));
    }
    ploidy = n;
}

bool Genotype::set(int v)
{
    std::pair< std::set<int>::iterator, bool > r = this->insert(v);
    if (ploidy < altSize())
    {
        Ploidy p(ploidy);
        std::stringstream message;
        message << (boost::format("*** Could not set allele number '%d' in %s event ***") % v % p.label()).str() << std::endl;
        message << "    The following alleles have already been used:" << std::endl;
        message << dump(message);
        BOOST_THROW_EXCEPTION(common::PreConditionException( message.str() ));
    }
    return r.second;
}

bool Genotype::reset(int v)
{
    return 0 < this->size() && 0 != this->erase(v);
}


std::ostream& operator<<( std::ostream& os, const Genotype& gt )
{
    //for (unsigned int i = gt.minPloidy(); i <= gt.maxPloidy(); i++)
    for (unsigned int i = 1;  // std::min(1,gt.minPloidy())
         i <= std::max(gt.getPloidy(),gt.maxPloidy());
         i++)
    {
        //if (i != gt.minPloidy())
        if (1 != i)  os << (gt.phased ? "|" : "/");
        os << int( gt.find(i) != gt.end() );
    }
    return os;
}

std::istream& operator>>( std::istream& is, Genotype& gt )
{
    std::string str;
    if (!(is >> str))
    {
        BOOST_THROW_EXCEPTION(common::PreConditionException("Could not parse GenoType information"));
    }
    std::vector<std::string> alleles;
    boost::split( alleles, str, boost::is_any_of("/|"));
    if (!alleles.empty())
        gt.setPloidy(alleles.size());
    for (unsigned int i = 0; i < alleles.size(); i++)
    {
        if (gt.altGtIndex_ == boost::lexical_cast<unsigned int>(alleles[i]))
            gt.set(i+1);
    }
    return is;
}


} // namespace model
} // namespace eagle
