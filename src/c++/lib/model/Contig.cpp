/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Mauricio Varea
 **/

#include <numeric>
#include <algorithm>

#include "model/Contig.hh"

namespace eagle
{
namespace model
{


void Contig::ins(Locus& loc, std::vector<char> seq)
{
    unsigned long pos = loc.pos();
    assert(id() == loc.chr());
    assert(size() >= pos);
    this->insert( begin()+pos, seq.begin(), seq.end() );
}


void Contig::del(Locus& loc, unsigned long len)
{
    unsigned long pos = loc.pos();
    assert(id() == loc.chr());
    assert(size() >= (pos + len));
    if (len)
    {
        this->erase( begin()+pos, begin()+pos+len );
    } else {
        this->erase( begin()+pos );
    }
}


unsigned long Contig::append(const std::vector<char>& seq, bool neg)
{
    std::vector<char> norm(seq.size());
    std::transform(seq.begin(), seq.end(), norm.begin(), IUPAC(neg));
    unsigned long size = this->size();
    this->insert( this->end(), norm.begin(), norm.end() );
    return static_cast<unsigned long>(this->size() - size);
}


// positions are: {1,2,3,... i,i+1,... j,j+1,...,-1}
// fetches [pos1,pos2[ bases if pos1<pos2, and "reversed ]pos2,pos1]" if pos1>pos2
// if pos2==-1, go until end
std::vector<char> Contig::read(const long pos1, const long pos2) const
{
    long p1(pos1);
    long p2(pos2);

    assert(p1 > 0);
    if (p2 == -1)
    {
        p2 = size()+1;
    }

    if (p1==p2)
    {
        return std::vector<char>();
    }
    else if (p1 < p2)
    {
        --p1;
        --p2;
        if ( (long)size() < p2 ) { std::cerr << "size<p2! (size=" << size() << ", p2=" << p2 << ")" << std::endl; assert(false); }
        assert( 0 <= p1 );
        return std::vector<char>( begin()+p1, begin()+p2 );
    }
    else if (p1 > p2)
    {
        if (p2 == -1)
        {
            p2 = 0;
        }
        assert( (long)size() >= p1 );
        assert( 0 <= p2 );
        return std::vector<char>( rbegin()+(size()-p1),
                                  rbegin()+(size()-p2) );
    }
    return std::vector<char>();
}


bool Contig::operator==( const Contig& other ) const
{
    return (this->id() == other.id());
}


} // namespace model
} // namespace eagle
