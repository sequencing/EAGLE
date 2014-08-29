/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Top level component that deals with variants.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_GENOME_VARIANT_LIST_HH
#define EAGLE_GENOME_VARIANT_LIST_HH

#include <string>
#include <vector>
#include <utility>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>

#include "io/Vcf.hh"
#include "genome/Event.hh"

namespace eagle
{
namespace genome
{


class VariantList: boost::noncopyable
{
public:
    VariantList(const std::vector<boost::filesystem::path> inputFiles,
                const boost::filesystem::path outputFile,
                const eagle::model::Ploidy ploidy,
                const bool overwrite = false);
    VariantList(const eagle::model::Ploidy ploidy);
    VariantList(const unsigned int ploidyLevel);

    void load( bool filterSnpsOut = false, bool filterBeginEndMarkersOut = false );
    unsigned int save(unsigned int i);
    void sort();
    void updateFirstEventPositionPerContig();
    boost::filesystem::path inputFile(unsigned int i) const  {return reader_.file(i);}
    boost::filesystem::path outputFile(unsigned int i) const {return writer_.file(i);}

    void pop_back()              {events_.pop_back();}
    void push_back(Event event)  {events_.push_back(event);}

    std::vector< Event >::const_iterator begin() const {return events_.begin();}
    std::vector< Event >::iterator begin() {return events_.begin();}
    std::vector< Event >::const_iterator end() const {return events_.end();}
    std::vector< Event >::iterator end() {return events_.end();}
    Event operator[](int i) const {return events_[i];}
    Event& operator[](int i) {return events_[i];}

    bool empty() const        {return events_.empty();}
    size_t size() const       {return events_.size();}
    size_t fileCount() const  {return reader_.size();}

    EventIterator findFirstEventForChromosome( std::string chr, eagle::model::Direction dir );
    void chromosomeNameCheck( const std::vector<std::string>& allContigNames );
    void removeDuplicatedTranslocations( );
    void pairing();
    bool pairingSearchInRange( const EventIterator& it, const EventIterator& first, const EventIterator& last );
    void check( const bool throwErrorIfTranslocationNotApplied );

    eagle::model::Ploidy ploidy() const {return ploidy_;}
protected:
    eagle::model::Ploidy ploidy_;
    eagle::io::VcfReader reader_;
    eagle::io::VcfWriter writer_;
    // events are stored in the canonical form (StructuralVariant, as opposed to VcfVariant)
    std::vector< Event > events_;
    std::vector< std::pair< std::string, EventIterator > > firstEventPositionPerContig_;
};


} // namespace genome
} // namespace eagle

#endif // EAGLE_GENOME_VARIANT_LIST_HH
