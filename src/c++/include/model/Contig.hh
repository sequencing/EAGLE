/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Definition of a contig
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_MODEL_CONTIG_HH
#define EAGLE_MODEL_CONTIG_HH

#include <cassert>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <utility>
#include <boost/filesystem.hpp>
#include "model/Struct.hh"
#include "model/Nucleotides.hh"


namespace eagle
{
namespace model
{


class Contig : public std::vector<char>
{
  public:
    Contig(std::string head) : name_(head) {}
    Contig() {}
    void reset() {name_=""; clear();}

    std::string name() const {return name_;}                                                                     // Used to be:
    std::string id() const {return( (name_.find_first_of(" ") != std::string::npos)                                 // .find_first_of(" |")
                                  ? std::string( name_.begin(), name_.begin() + name_.find_first_of(" ") )          // .find_first_of(" |")
                                  : name_);}                        // returns entire name if delim not found
    std::string remainder() const {return( (name_.find_first_of(" ") != std::string::npos)                          // .find_first_of(" |")
                                         ? std::string(name_.begin() + name_.find_first_of(" ") + 1, name_.end() )  // .find_first_of(" |")
                                         : std::string() );}        // returns empty if delim not found
    void name(std::string head) {name_ = head;}
    void name(std::string id,std::string remainder) {name_ = (remainder.length() ? (id + " " + remainder) : id );}  // (id + "|" + remainder)

    char get(Locus& loc) const
    { assert(id() == loc.chr()); return *(begin() + loc.pos() - 1); }
    void put(const char base, const bool neg = false)
    { std::vector<char>::push_back( iupac.normalize(base,neg) ); }

    void ins(Locus& loc, std::vector<char> seq);
    void del(Locus& loc, unsigned long len = 0);
    unsigned long append(const std::vector<char>& seq, const bool neg = false);

    // read: positions are: {1,2,3,... i,i+1,... j,j+1,...,-1}
    // fetches [pos1,pos2[ bases if pos1<pos2, and "reversed ]pos2,pos1]" if pos1>pos2
    // if pos2==-1, go until end
    std::vector<char> read(const long pos1=1, const long pos2=-1) const;

    bool operator==( const Contig& other ) const;

  private:
    std::string name_;
    IUPAC iupac;
};



} // namespace model
} // namespace eagle

#endif // EAGLE_MODEL_CONTIG_HH
