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

#ifndef EAGLE_MODEL_GENOTYPE_HH
#define EAGLE_MODEL_GENOTYPE_HH

#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

namespace eagle
{
namespace model
{


class Ploidy
{
public:
    Ploidy( unsigned int level, std::map<std::string,unsigned int> exceptions )
    : organismPloidy_(level), nonPloidyChromosomes_(exceptions) {}
    Ploidy( unsigned int level = 0)
    : organismPloidy_(level) {}

    unsigned int level(const std::string& chrom) const;
    unsigned int max() const;
    std::string label(unsigned int p) const;
    std::string label() const {return label(organismPloidy_);}
private:
    unsigned int organismPloidy_;
    std::map<std::string,unsigned int> nonPloidyChromosomes_;
};

/*
class Genotype : public std::vector<bool>
{
public:
    // initialize Genotype as Homo-ref
    Genotype(unsigned int n)
    : std::vector<bool>(n,false), phased_(false) {}
    // assume haploid, if no ploidy is given
    Genotype()
    : std::vector<bool>(1,false), phased_(false) {}

    // an event cannot be Homo-ref, since this implies a variant has not been applied
    bool isHomozygousRef() const {return this->size() == (size_t) std::count(this->begin(),this->end(),false);}
    void setHomozygousRef() {std::fill(this->begin(),this->end(),false);}  // not processed
    // non-conflicting events are Homo-diff by default
    bool isHomozygousDiff() const {return this->size() == (size_t) std::count(this->begin(),this->end(),true);}
    void setHomozygousDiff() {std::fill(this->begin(),this->end(),true);}  // processed in all alleles
    // any of the above
    bool isHomozygous() const {return isHomozygousDiff() || isHomozygousRef();}
    // conflicting events end up as Het
    bool isHeterozygous() const {return !isHomozygous();}

public:
    bool phased_;
};
*/


/*
 * \brief Genotype contains discrete set of alleles
 *    -1: any allele
 *     0: REF allele
 *  1..n: ALT alleles
 */
class Genotype : public std::set<int>
{
public:
    Genotype(unsigned int n = 1, unsigned int altGtIndex = 1)
        : phased(false), ploidy(n), altGtIndex_(altGtIndex) {}

    void setPloidy(unsigned int n);
    unsigned int getPloidy() const {return ploidy;}

    std::set<int>::iterator altBegin() const {int i=1; while (this->find(--i) == this->end()) {if (-2 == i) return this->begin();} return this->upper_bound(i);}
    std::set<int>::iterator altBegin() {int i=1; while (this->find(--i) == this->end()) {if (-2 == i) return this->begin();} return this->upper_bound(i);}

    // a processed event cannot be Homo-ref, since this implies a variant has not been applied
    //bool isHomozygousRef() const {return this->empty() || this->rbegin().base() == this->find(0);}
    bool isHomozygousRef() const  {return this->empty() || ( 1 == this->size() && 0 == *(this->begin()) );}
    bool isHomozygousDiff() const {return             ploidy == this->size() && 0 < *(this->begin());}
    bool isHeterozygous() const   {return                  1 == this->size() && 0 < *(this->begin());}
    bool isHomozygous() const {return isHomozygousDiff() || isHomozygousRef();}

    unsigned int minPloidy() const {return (0 == this->size() || 0 >= *(this->rbegin())) ? 1      : static_cast<unsigned int>( *(altBegin()) );}
    unsigned int maxPloidy() const {return (0 == this->size() || 0 >= *(this->rbegin())) ? ploidy : static_cast<unsigned int>( *(this->rbegin()) );}

    unsigned int altSize() const {return this->size() - static_cast<size_t>(this->find(0) != this->end())
                                                      - static_cast<size_t>(-1 == *(this->begin()));}

    bool set(int v);
    bool reset(int v);
    std::ostream& dump(std::ostream& os);
    friend std::istream& operator>>( std::istream& is, Genotype& gt );

    bool phased;

private:
    unsigned int ploidy;
    unsigned int altGtIndex_;
};


std::ostream& operator<<( std::ostream& os, const Genotype& gt );


} // namespace model
} // namespace eagle

#endif // EAGLE_MODEL_GENOTYPE_HH
