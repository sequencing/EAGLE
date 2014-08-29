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

#ifndef EAGLE_MODEL_STRUCT_HH
#define EAGLE_MODEL_STRUCT_HH

#include "common/Logger.hh"
#include "model/Nucleotides.hh"

#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


namespace eagle
{
namespace model
{


struct Direction
{
    enum DirectionType { NONE=0, FWD=1, REV=2, BIDIR=3 } value_;
    Direction( enum DirectionType val )
    {
        value_ = val;
    }

    int offset(bool defaultIsFwd=false) const
    {
        switch (value_)
        {
        case FWD: return 1;
        case REV: return -1;
        case BIDIR: return static_cast<int>(defaultIsFwd);
        default:
            if (defaultIsFwd) return 1;
            else return 0; //EAGLE_ERROR("Unexpected direction");
        }
    }

    std::string str() const
    {
        switch (value_)
        {
        case NONE:  return "";
        case FWD:   return ">";
        case REV:   return "<";
        case BIDIR: return "<>";
        }
        EAGLE_ERROR("Unexpected direction");
    }

    bool defined() const { return (NONE  != value_); }
    bool isFwd() const   { return (FWD   == value_); }
    bool isRev() const   { return (REV   == value_); }
    bool isBiDir() const { return (BIDIR == value_); }
    bool operator==( const Direction& rhs ) const { return (value_ == rhs.value_); }
    bool operator!=( const Direction& rhs ) const { return (value_ != rhs.value_); }
    bool operator<( const Direction& rhs ) const  { return (value_ < rhs.value_); }

    bool sameAs( const Direction& rhs ) const
    {
        switch (value_)
        {
        case NONE:
        case BIDIR:
             return true;
        case FWD:   return (rhs.value_ != REV);
        case REV:   return (rhs.value_ != FWD);
        }
        EAGLE_ERROR("Unexpected direction");
    }

    enum DirectionType inv() const
    {
        switch(value_)
        {
        case NONE:  return BIDIR;
        case FWD:   return REV;
        case REV:   return FWD;
        case BIDIR: return NONE;
        }
        EAGLE_ERROR("Unexpected direction");
    }

};


class Locus
{
public:
    Locus(std::string s,unsigned long u, bool half=false) : chr_(s), pos_(2 * u + static_cast<unsigned long>(half)) {}
    Locus(std::string obj);
    Locus() {};
    virtual ~Locus() {}

    static bool lt( const Locus& lhs, const Locus& rhs );
    static bool eq( const Locus& lhs, const Locus& rhs );
    static bool gt( const Locus& lhs, const Locus& rhs );
    bool operator<( const Locus & rhs ) const;
    bool operator>( const Locus & rhs ) const;
    bool operator==( const Locus& rhs ) const { return chr_ == rhs.chr() && pos_ == static_cast<unsigned long>(2 * rhs.decimalPos()); }
    bool operator!=( const Locus& rhs ) const { return !(*this == rhs); }
    Locus& operator+=(const Locus& rhs) { pos_ += static_cast<unsigned long>(2 * rhs.decimalPos()); return *this; }
    Locus& operator-=(const Locus& rhs) { assert(pos_ >= static_cast<unsigned long>(2 * rhs.decimalPos()) && "Larger delta than allowed would result in negative positions");
                                          pos_ -= static_cast<unsigned long>(2 * rhs.decimalPos()); return *this; }

    void pos(unsigned long p, bool half=false) { pos_ = 2 * p + static_cast<unsigned long>(half); }
    virtual unsigned long pos() const { return pos_ / 2; }
    double decimalPos() const { return static_cast<double>(pos_) / 2; }
    void chr(std::string c) { chr_ = c; }
    std::string chr() const { return chr_; }

protected:
    std::string chr_;
    unsigned long pos_;    // pos_ = 2 * real_pos_in_genome
};

std::ostream& operator<<( std::ostream& os, const Locus& loc );


class Breakend : public Locus
{
  public:
    Breakend(std::string s, unsigned long u, Direction d=Direction::NONE, std::string b=".") : Locus(s,u), dir(d), base(b) {}
    Breakend(std::string obj, Direction d=Direction::NONE, std::string b=".") : Locus(obj), dir(d), base(b) {}
    virtual ~Breakend() {}

    std::string direction() const;
    bool hasSameLocus( const Breakend& rhs ) const { return Locus::operator==(dynamic_cast<const Locus&>(rhs)); }
    bool lessThanLocusComparison( const Breakend& rhs ) const { return Locus::operator<(dynamic_cast<const Locus&>(rhs)); }
    bool operator==( const Breakend& rhs ) const { return this->hasSameLocus(rhs)
                                                       && this->dir == rhs.dir; }
    void pos(unsigned long p, bool half=false) {Locus::pos(p,half);}
    unsigned long pos(Direction d) const { return d.isRev()
                                                ? (unsigned long)ceil(decimalPos())
                                                : (unsigned long)floor(decimalPos());
    }
    unsigned long posAfter(Direction d) const { return d.isRev()
                                                ? (unsigned long)ceil(decimalPos()-1)
                                                : (unsigned long)floor(decimalPos()+1);
    }
    virtual unsigned long pos() const { return pos(dir); }
    virtual unsigned long posAfter() const { return posAfter(dir); }
    bool operator<( const Breakend & rhs ) const
    {
        if (Locus::operator<(dynamic_cast<const Locus&>(rhs)))
            return true;
        if (Locus::operator>(dynamic_cast<const Locus&>(rhs)))
            return false;
        return (dir < rhs.dir);
    }
    Breakend& operator+=(const Locus& rhs)
    {
        if (dir.isRev())
        {
            return dynamic_cast<Breakend&>(Locus::operator-=(rhs));
        } else {
            return dynamic_cast<Breakend&>(Locus::operator+=(rhs));
        }
    }
    Breakend& operator-=(const Locus& rhs)
    {
        if (dir.isRev())
        {
            return dynamic_cast<Breakend&>(Locus::operator+=(rhs));
        } else {
            return dynamic_cast<Breakend&>(Locus::operator-=(rhs));
        }
    }
    Direction dir;
    std::string base;
};

std::ostream& operator<<( std::ostream& os, const Breakend& bnd );
// shift-right for output streams: Very Nasty!
//std::ostream& operator>>( std::ostream& os, const Breakend& bnd );

/*
 * \brief Elementary unit of a Structural Variant
 *
 * Strictly speaking:
 *     Adjacency = pair<Breakend,Breakend>; Mates = pair<Adjacency,Adjacency>; ComplexRearrangement = Mates+
 * But we do:
 *     Adjacency = pair<Breakend,Breakend>; ComplexRearrangement = Adjacency+
 * And later pair the <Event>s
 */
struct ComplexRearrangement
{
    ComplexRearrangement(Breakend bnd1, Breakend bnd2, std::string seq="", unsigned int altGtIndex=1)
    : adjacency(std::make_pair(bnd1,bnd2))
    , sequence(seq.begin(),seq.end())
    , altGtIndex_(altGtIndex)
    {
        if (bnd1.dir.isRev() && !sequence.empty())
        {
            std::transform(sequence.begin(), sequence.end(), sequence.begin(), IUPAC(true));
            std::reverse(sequence.begin(), sequence.end());
        }
    }

    void inverse();
    void setDirection(Direction d1,Direction d2) {this->adjacency.first.dir = d1; this->adjacency.second.dir = d2;}
    void setDirection(Direction d) {this->setDirection(d,d);}
    void setBase(std::string b1,std::string b2) {this->adjacency.first.base = b1; this->adjacency.second.base = b2;}
    void setBase(std::string b) {this->setBase(b,b);}
    std::pair<Breakend,Breakend> adjacency;
    std::vector<char> sequence;
    unsigned int altGtIndex_;
};

std::ostream& operator<<( std::ostream& os, const ComplexRearrangement& obj );
// shift-right for output streams: Very Nasty!
//std::ostream& operator>>( std::ostream& os, const ComplexRearrangement& obj );



} // namespace model
} // namespace eagle

#endif // EAGLE_MODEL_STRUCT_HH
