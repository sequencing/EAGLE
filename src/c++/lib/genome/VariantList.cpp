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

#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "model/Nucleotides.hh"
#include "genome/VariantList.hh"

using namespace std;
using namespace boost::assign;
using namespace boost::lambda;


namespace eagle
{
namespace genome
{


VariantList::VariantList (const vector<boost::filesystem::path> inputFiles,
                          const boost::filesystem::path outputFile,
                          const eagle::model::Ploidy ploidy,
                          const bool overwrite)
    : ploidy_( ploidy )
    , reader_( inputFiles )
    , writer_( outputFile, overwrite )
    , events_()
{
    for (vector<boost::filesystem::path>::const_iterator var = inputFiles.begin(); var != inputFiles.end(); ++var)
    {
        clog << "+ Input variant list: " << *var << endl;
    }
    for (vector<boost::filesystem::path>::const_iterator var = writer_.begin(); var != writer_.end(); ++var)
    {
        clog << "+ Output variant list: " << *var << endl;
    }
}
VariantList::VariantList (const eagle::model::Ploidy ploidy)
    : ploidy_( ploidy )
    , events_()
{}
VariantList::VariantList (const unsigned int ploidyLevel)
    : ploidy_( ploidyLevel )
    , events_()
{}


void VariantList::load( bool filterSnpsOut, bool filterBeginEndMarkersOut )
{
    eagle::io::VcfVariant variant;
    while (reader_.getNextVariant(variant, filterSnpsOut, filterBeginEndMarkersOut))
    {
        for(eagle::io::VcfVariant::const_iterator v=variant.begin(); v!=variant.end(); ++v)
        {
            if (filterSnpsOut && v->getType() == eagle::model::variant::SNP)
            {
                continue;
            }
            if (filterBeginEndMarkersOut && v->getType() == eagle::model::variant::Undefined)
            {
                continue;
            }
            unsigned int thisPloidy = ploidy_.level(v->getVariant().adjacency.first.chr());
            EAGLE_DEBUG( 0, "... " << *v ); // << " (P = " << thisPloidy << ")");
            if (v->getVariant().adjacency.first.dir == eagle::model::Direction::NONE
            ||  v->getVariant().adjacency.second.dir == eagle::model::Direction::NONE)
            {   // create forward variant
                eagle::model::StructuralVariant fwd(*v);

                if (eagle::model::variant::SNP == v->getType())
                {
                    fwd.getVariant().setDirection(eagle::model::Direction::BIDIR);
                }
                else
                {
                    fwd.getVariant().setDirection(eagle::model::Direction::FWD);
                }

                EAGLE_DEBUG( 0, "..... " << fwd );
                events_.push_back(Event( fwd, variant.metadata_, thisPloidy ));

                if (v->hasDeletion() || v->hasInsertion())
                {
                    // create reverse variant
                    eagle::model::StructuralVariant rev(*v);
                    rev.getVariant().inverse();
                    EAGLE_DEBUG( 0, "..... " << rev );
                    events_.push_back(Event( rev, variant.metadata_, thisPloidy ));
                }
            } else { // translocation
                EAGLE_DEBUG( 0, "..... " << *v );
                events_.push_back(Event( *v, variant.metadata_, thisPloidy ));
            }
        }
    }
}

unsigned int VariantList::save(unsigned int i)
{
    EventIterator event(events_.begin());
    unsigned int eventCount(0);
    writer_.open(i);
    writer_.writeHeader();
    while(events_.end() != event)
    {
        // we don't want to write the begin/end markers
        // we do - if (event->isDefined())
        {
            if (!event->allele_.isHomozygousRef())
            {   // we only want to show in screen the ones that have actually been processed (file may contain more)
                EAGLE_DEBUG( 0, "... [processed] " << *event );
            }
            writer_.write(eagle::io::VcfVariant(event->getStructuralVariant(),event->metadata_));
            ++eventCount;
        }
        ++event;
    }
    return eventCount;
}

void VariantList::sort()
{
    stable_sort( events_.begin(), events_.end() );
    updateFirstEventPositionPerContig();
}

void VariantList::updateFirstEventPositionPerContig()
{
    // update index per chromosome name
    firstEventPositionPerContig_.clear();
    string lastChrName( "a contig name that will never happen" );
    for (EventIterator it = events_.begin(); it != events_.end(); ++it)
    {
        // Test source chromosome name if it differs from the last one tested
        if (it->src() != lastChrName)
        {
            lastChrName = it->src();
            firstEventPositionPerContig_.push_back( make_pair( lastChrName, it ) );
        }
    }
}

EventIterator VariantList::findFirstEventForChromosome( string chr, eagle::model::Direction dir )
{
    for( vector< pair< string, EventIterator > >::iterator it = firstEventPositionPerContig_.begin() ;
         it != firstEventPositionPerContig_.end() ;
         ++it)
    {
        if ( it->first == chr )
        {
            if ( dir.isFwd() )
            {
                return it->second;
            }
            else if ( dir.isRev() )
            {
                ++it;
                if (it != firstEventPositionPerContig_.end())
                {
                    return it->second - 1;
                }
                else
                {
                    return events_.end() - 1;
                }
            }
            else {
                EAGLE_ERROR("Unexpected direction");
            }
        }
    }
    return events_.end();
}

void VariantList::chromosomeNameCheck( const vector<string>& allContigNames )
{
    string lastChrName( "a contig name that will never happen" );
    for (EventIterator it = events_.begin(); it != events_.end(); ++it)
    {
        // Test source chromosome name if it differs from the last one tested
        if (it->src() != lastChrName)
        {
            EAGLE_DEBUG(5, "[checking] looking for chromosome name: " << it->src() );
            if (find( allContigNames.begin(), allContigNames.end(), it->src()) == allContigNames.end())
            {
                EAGLE_ERROR( (boost::format("Unexpected chromosome name in variants file: %s") % it->src()).str() );
            }

            lastChrName = it->src();
        }

        // Test dest chromosome name if it differs from the source
        // (This wastes a bit of time, but only happens for translocations, which shouldn't be too frequent)
        if (lastChrName != it->dest())
        {
            EAGLE_DEBUG(5, "[checking] looking for chromosome name: " << it->dest() );
            if (find( allContigNames.begin(), allContigNames.end(), it->dest()) == allContigNames.end())
            {
                EAGLE_ERROR( (boost::format("Unexpected chromosome name in variants file: %s") % it->dest()).str() );
            }
        }
    }
}

void VariantList::removeDuplicatedTranslocations( )
{
    EventIterator it = events_.begin();
    while (it != events_.end() && (it+1) != events_.end())
    {
        EventIterator nextIt = it+1;

        if (it->hasTranslocation() && nextIt->hasTranslocation() && it->src() != it->dest() && *it == *nextIt )
        {
            EAGLE_WARNING( "Removing duplicated translocation: " << *it );
            events_.erase( nextIt );
        }
        else
        {
            ++it;
        }
    }

    updateFirstEventPositionPerContig();
}

static bool ltComparePotentialPairByPosition( const Event& event1, const Event& event2 )
{
    EAGLE_DEBUG(5, "[pairing test] event1: " << event1.getStructuralVariant() );
    EAGLE_DEBUG(5, "[pairing test] event2: " << event2.getStructuralVariant() );
    bool result = event1.getVariant().adjacency.first.lessThanLocusComparison( event2.getVariant().adjacency.first );
    EAGLE_DEBUG(5, "[pairing test] result: " << result );
    return result;
}

void VariantList::pairing()
{
    for (EventIterator it = events_.begin(); it != events_.end(); ++it)
    {
        if (!it->incoming().defined()) { continue; } // skip chromosome extremities 'fake' events
        if (it->incoming().isBiDir())
        {
            EAGLE_DEBUG(5, "[pairing] skipping the following (bi-directional) event: " << it->getStructuralVariant() );
            it->pairedEvent_ = it - events_.begin(); // Pair bi-directional events with themselves
            continue;
        }
        if (it->pairedEvent_ == 0 && it != events_.end()-1)
        {
            EAGLE_DEBUG(5, "[pairing] looking for paired event for: " << it->getStructuralVariant() );
            bool foundPairedEvent = false;

            // Linear local search
            {
                const unsigned int maxLinearSearchDistance = 5;
                EventIterator first = it+1;
                EventIterator last  = events_.end();
                if (last - first > maxLinearSearchDistance)
                {
                    last = first + maxLinearSearchDistance;
                }
                foundPairedEvent = pairingSearchInRange( it, first, last );
            }

            // Global search
            if (!foundPairedEvent)
            {
                Event fakeComplementaryEvent = *it;
                fakeComplementaryEvent.getVariant().adjacency.first = it->getVariant().adjacency.second;
                pair< EventIterator, EventIterator > boundaries = equal_range( it+1, events_.end(), fakeComplementaryEvent, ltComparePotentialPairByPosition );
                foundPairedEvent = pairingSearchInRange( it, boundaries.first, boundaries.second );
            }

            if (!foundPairedEvent)
            {
                stringstream message;
                message << "*** Could not pair the following event ***" << endl
                        << "    " << *it << endl;
                BOOST_THROW_EXCEPTION(common::PreConditionException( message.str() ));
            }
        }
        else
        {
            // EAGLE_DEBUG(0, "..... already known paired event for " << *it);
        }
    }
}

bool VariantList::pairingSearchInRange( const EventIterator& it, const EventIterator& first, const EventIterator& last )
{
    bool foundPairedEvent = false;
    for (EventIterator it2 = first; it2 != last; ++it2)
    {
        EAGLE_DEBUG(5,"[pairing] ... trying event: " << it2->getStructuralVariant() );
        if (it2->pairedEvent_) { continue; }
        if (it->getStructuralVariant().getType() != it2->getStructuralVariant().getType()) { continue; }

        if ( (it2->getVariant().adjacency.first.hasSameLocus( it->getVariant().adjacency.second ))
             && (it2->getVariant().adjacency.second.hasSameLocus( it->getVariant().adjacency.first ))
            )
        {
            if (it->getStructuralVariant().hasSNP() && (it->incoming().isFwd() || it->incoming().isRev()) )
            {
                foundPairedEvent = true;
            } else if ( (it2->incoming() != it->outgoing()) ) {
                foundPairedEvent = true;
            }
        }
        if (foundPairedEvent)
        {
            EAGLE_DEBUG(5,"[pairing] ... found paired event: " << it2->getStructuralVariant() );
            it->pairedEvent_ = it2 - events_.begin();
            it2->pairedEvent_ = it - events_.begin();
            EAGLE_DEBUG(5,"[pairing]           at distance: " << it2-it );
            break;
        }
    }
    return foundPairedEvent;
}

void VariantList::check( const bool throwErrorIfTranslocationNotApplied )
{
    clog << "Checking if some events were not applied..." << endl;
    int notAppliedTranslocationCount = 0;
    int appliedTranslocationCount = 0;
    for (EventIterator it = events_.begin(); it != events_.end(); ++it)
    {
        if (!it->incoming().defined()) { continue; } // skip fake events
        if (it->allele_.isHomozygousRef())
        {
            clog << "... event not applied: " << *it << endl;
            it->allele_.set(-1);  // mark as processed, in no particular allele

            if (it->pairedEvent_)
            {
                clog << "..... paired event: " <<  events_[it->pairedEvent_] << endl;
                events_[it->pairedEvent_].allele_.set(-1);  // mark as processed, in no particular allele
            }
            if (it->hasTranslocation())
            {
                ++notAppliedTranslocationCount;
                EAGLE_WARNING("The above translocation has not been applied");
            }
        }
        else if (it->hasTranslocation())
        {
            ++appliedTranslocationCount;
        }
    }
    if (notAppliedTranslocationCount)
    {
        if (throwErrorIfTranslocationNotApplied)
        {
            EAGLE_ERROR( (boost::format("%d translocation(s) were not applied") % notAppliedTranslocationCount).str() );
        }
        else
        {
            EAGLE_WARNING( (boost::format("%d translocation(s) were not applied") % notAppliedTranslocationCount).str() );
        }
    }
    else if (appliedTranslocationCount)
    {
        clog << "\tAll translocations applied!" << endl;
    }
}

} // namespace genome
} // namespace eagle
