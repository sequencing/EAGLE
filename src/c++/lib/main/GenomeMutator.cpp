/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Top level component to induce variants in a reference.
 **
 ** \author Mauricio Varea
 **/

#include <iostream>
#include <utility>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#ifdef DISTRIBUTED_GENOME_MUTATOR
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include "genome/SharedFastaReference.hh"
#endif //ifdef DISTRIBUTED_GENOME_MUTATOR

#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "genome/SharedFastaReference.hh"
#include "model/Struct.hh"
#include "model/Genotype.hh"
#include "main/GenomeMutator.hh"


using namespace std;


namespace eagle
{
namespace main
{


#ifndef DISTRIBUTED_GENOME_MUTATOR

void GenomeMutator::run()
{
    if (options_.onlyPrintOutputContigNames)
    {
        BOOST_FOREACH( const string& contigName, genome_.allContigNames() )
        {
            unsigned int ploidy = variantList_.ploidy().level(contigName);
            for (unsigned int i=0; i<ploidy; ++i)
            {
                string outName = prefixToAdd_ + contigName
                    + (boost::format("_allele%d") % (i+1)).str()
//                    + ((direction.isFwd()) ? ""                                         : "_rev")
                    ;
                cout << "Chromosome allele: " << outName << endl;
            }
        }
        return;
    }

    unsigned long timeProcessing = 0;
    unsigned long timeIO = 0;

    clock_t start = clock();

    std::clog << "Loading " << variantList_.fileCount() << " variant list" << (variantList_.fileCount() - 1 ? "s " : " ") << "..." << std::endl;
    variantList_.load();
    std::clog << "Loaded " << variantList_.size() << " event" << (variantList_.size() - 1 ? "s " : " ")
              << "in " << eagle::common::displayTime(clock() - start, timeIO) << std::endl;

    start = clock();
    std::clog << "Loading " << genome_.fileCount() << " reference genome" << (genome_.fileCount() - 1 ? "s " : " ") << "..." << std::endl;
    genome_.load();
    std::clog << "Loaded " << genome_.contigCount() << " contig" << (genome_.contigCount() - 1 ? "s " : " ")
              << "in " << eagle::common::displayTime(clock() - start, timeIO) << std::endl;

    std::clog << "Total genome size is " << genome_.length() << std::endl;

    std::map< eagle::genome::ReferenceIterator, unsigned int > forwardStartingPoints;
    std::map< eagle::genome::ReferenceIterator, unsigned int > reverseStartingPoints;
    // Adding chromosome begin&end markers
    for (eagle::genome::ReferenceIterator contig=genome_.begin(); contig != genome_.end(); ++contig)
    {
        unsigned int n( variantList_.ploidy().level(contig->id()) );
        forwardStartingPoints[contig] = n;
        reverseStartingPoints[contig] = n;
        eagle::model::StructuralVariant sv1( contig->id(), 0 );
        eagle::model::StructuralVariant sv2( contig->id(), contig->size()+1 );
        sv1.getVariant().adjacency.second.pos( sv1.getVariant().adjacency.first.pos() );
        sv2.getVariant().adjacency.second.pos( sv2.getVariant().adjacency.first.pos() );
        variantList_.push_back( eagle::genome::Event( sv1, n ));
        variantList_.push_back( eagle::genome::Event( sv2, n ));
    }

    start = clock();
    std::clog << "Sorting variant list..." << std::endl;
    variantList_.sort();
    //events_.resize( std::unique(variantList_.begin(),variantList_.end()) - variantList_.begin() );  // remove duplicates
    for (eagle::genome::EventIterator it = variantList_.begin(); it != variantList_.end(); ++it)
    {   // compiler should optimize this away in Release mode
        EAGLE_DEBUG( 0, "... " << it->getStructuralVariant() );
    }
    std::clog << "Sorted " << variantList_.size() << " event" << (variantList_.size() - 1 ? "s " : " ")
              << "in " << eagle::common::displayTime(clock() - start, timeProcessing) << std::endl;

    variantList_.chromosomeNameCheck( genome_.allContigNames() );

    // Check if some translocations appear multiple times (for example if the user specifies the same input file multiple times).
    // Warn and remove duplicated translocations
    variantList_.removeDuplicatedTranslocations();

    start = clock();
    std::clog << "Pairing variant list..." << std::endl;
    variantList_.pairing();
    std::clog << "Paired " << variantList_.size() << " event" << (variantList_.size() - 1 ? "s " : " ")
              << "in " << eagle::common::displayTime(clock() - start, timeProcessing) << std::endl;

    assert( forwardStartingPoints.size() > 0 && "No chromosome found");
    eagle::model::Direction direction = eagle::model::Direction::FWD;
    std::map< eagle::genome::ReferenceIterator, unsigned int >::iterator startingPoint = forwardStartingPoints.begin();
    while (startingPoint != reverseStartingPoints.end())
    {
        eagle::genome::ReferenceIterator contig = startingPoint->first;
        unsigned int n = startingPoint->second;

        if (n > 0)
        {
            clog << "Started processing " << n << " allele(s) from " << ((direction.isFwd())?"beginning":"end") << " of chromosome " << contig->id() << endl;

            sample_.clear();
            sample_.resize(n);
            for(eagle::genome::iterator contigOut=sample_.begin(); contigOut != sample_.end(); ++contigOut)
            {
                unsigned int i = static_cast<unsigned int>(contigOut - sample_.begin());
                contigOut->name( prefixToAdd_ + contig->id()
                                 + (boost::format("_allele%d") % (i+1)).str()
                                 + ((direction.isFwd()) ? ""                                         : "_rev")
                                 , contig->remainder() );

                eagle::model::Locus finalPosition(contig->id(),(direction.isFwd())?contig->size():0);
                start = clock();
                unsigned int eventCount;
                std::clog << "Applying structural variants to '" << contig->name() << "' allele " << (i+1) << "..." << std::endl;
                //                EventIterator startPosition( (direction.isRev())?(bounds.second-1):bounds.first );
                eagle::genome::EventIterator startPosition = variantList_.findFirstEventForChromosome( contig->id(), direction );
                assert( startPosition != variantList_.end() );

                eventCount = process(i+1, startPosition, *contigOut, direction, finalPosition );

                std::clog << "Applied " << eventCount << " event" << (eventCount==1 ? " " : "s ")
                          << "to a " << variantList_.ploidy().label(n) << " chromosome "
                          << "in " << eagle::common::displayTime(clock() - start, timeProcessing) << std::endl;

                // Remove the starting point that we just processed
                --(startingPoint->second);
                // Remove the final position that we just reached from lists of startingPoints
                std::map< eagle::genome::ReferenceIterator, unsigned int >& mapToSearch = (0 == finalPosition.pos())?forwardStartingPoints:reverseStartingPoints;
                for (std::map< eagle::genome::ReferenceIterator, unsigned int >::iterator it = mapToSearch.begin(); it != mapToSearch.end(); ++it) {
                    if (it->first->id() == finalPosition.chr()) {
                        EAGLE_DEBUG(5, "[starting points] Deleting starting point at " << ((0 == finalPosition.pos())?"beginning":"end") << " of chromosome " << finalPosition.chr());
                        if ( it->second <= 0 ) {
                            EAGLE_ERROR("[starting points] Trying to delete already-deleted starting points at " + std::string((1 == finalPosition.pos())?"beginning":"end") + " of chromosome " + finalPosition.chr());
                        }
                        --(it->second);
                        break;
                    }
                }
            }
            start = clock();
            std::clog << "Saving " << sample_.fileCount() << " sample chromosome" << (sample_.fileCount() > 1 ? "s " : " ") << "..." << std::endl;
            sample_.save();
            std::clog << "Saved " <<  sample_.contigCount() << " contig" << (sample_.contigCount() > 1 ? "s " : " ")
                      << "in " << eagle::common::displayTime(clock() - start, timeIO) << std::endl;
        }

        // Second loop will process things backwards
        ++startingPoint;
        if (startingPoint == forwardStartingPoints.end())
        {
            startingPoint = reverseStartingPoints.begin();
            direction = eagle::model::Direction::REV;
        }
    }
    start = clock();
    sample_.saveMetadata();
    std::clog << "Saved " GENOMESIZE_XML
              << " in " << eagle::common::displayTime(clock() - start, timeIO) << std::endl;

    start = clock();
    unsigned int totalVariantCount(0);
    // for (unsigned int i=0; i<variantList_.ploidy().max(); i++)
    for (unsigned int i=0; i<1; i++)
    {
        std::clog << "Saving " << variantList_.outputFile(i) << "..." << std::endl;
        unsigned int variantCount = variantList_.save(i);
        //std::clog << "..saved " << variantCount << " event" << (variantCount - 1 ? "s " : " ") << std::endl;
        totalVariantCount += variantCount;
    }
    std::clog << "Saved canonical vcf (" << totalVariantCount << " event" << (totalVariantCount - 1 ? "s" : "")
              << ") in " << eagle::common::displayTime(clock() - start, timeIO) << std::endl;

    std::clog << "+ Total processing time: " << eagle::common::displayTime(timeProcessing) << std::endl;
    std::clog << "+ Total I/O time: " << eagle::common::displayTime(timeIO) << std::endl;
    variantList_.check( !options_.noTranslocationError );
}

#else //ifndef DISTRIBUTED_GENOME_MUTATOR

void GenomeMutator::run()
{
    unsigned long timeProcessing = 0;
    unsigned long timeIO = 0;

    clock_t start = clock();

    clog << "Loading " << variantList_.fileCount() << " variant list" << (variantList_.fileCount() - 1 ? "s " : " ") << "..." << endl;
    const bool filterSnpOut = true;
    variantList_.load( filterSnpOut );
    clog << "Loaded " << variantList_.size() << " event" << (variantList_.size() - 1 ? "s " : " ")
              << "in " << eagle::common::displayTime(clock() - start, timeIO) << endl;

    start = clock();
    clog << "Loading " << genome_.fileCount() << " reference genome" << (genome_.fileCount() - 1 ? "s " : " ") << "..." << endl;
//    genome_.load();
    const vector<string> contigIds            = eagle::genome::SharedFastaReference::get()->allContigNames();
    const vector<unsigned long> contigLengths = eagle::genome::SharedFastaReference::get()->allContigLengths();
    assert( contigIds.size()   == contigLengths.size() );
    clog << "Found " << contigIds.size() << " contig(s)." << endl;
    assert( contigIds.size() > 0 && "No chromosome found");
    clog << "Loaded " << genome_.contigCount() << " contig" << (genome_.contigCount() - 1 ? "s " : " ")
              << "in " << eagle::common::displayTime(clock() - start, timeIO) << endl;

//    clog << "Total genome size is " << genome_.length() << endl;

    vector<unsigned int > forwardStartingPoints( contigIds.size() );;
    vector<unsigned int > reverseStartingPoints( contigIds.size() );;
    // Adding chromosome begin&end markers
/*
    for (eagle::genome::ReferenceIterator contig=genome_.begin(); contig != genome_.end(); ++contig)
    {
        unsigned int n( variantList_.ploidy().level(contig->id()) );
        forwardStartingPoints[contig] = n;
        reverseStartingPoints[contig] = n;
        eagle::model::StructuralVariant sv1( contig->id(), 0 );
        eagle::model::StructuralVariant sv2( contig->id(), contig->size()+1 );
        sv1.getVariant().adjacency.second.pos( sv1.getVariant().adjacency.first.pos() );
        sv2.getVariant().adjacency.second.pos( sv2.getVariant().adjacency.first.pos() );
        variantList_.push_back( eagle::genome::Event( sv1, n ));
        variantList_.push_back( eagle::genome::Event( sv2, n ));
    }
*/
    for (unsigned int i=0; i<contigIds.size(); ++i)
    {
        const string&       contigId     = contigIds[i];
        const unsigned long contigLength = contigLengths[i];

        unsigned int n( variantList_.ploidy().level(contigId) );
        forwardStartingPoints[i] = n;
        reverseStartingPoints[i] = n;
        eagle::model::StructuralVariant sv1( contigId, 0 );
        eagle::model::StructuralVariant sv2( contigId, contigLength+1 );
        sv1.getVariant().adjacency.second.pos( sv1.getVariant().adjacency.first.pos() );
        sv2.getVariant().adjacency.second.pos( sv2.getVariant().adjacency.first.pos() );
        variantList_.push_back( eagle::genome::Event( sv1, n ));
        variantList_.push_back( eagle::genome::Event( sv2, n ));
    }



    start = clock();
    clog << "Sorting variant list..." << endl;
    variantList_.sort();
    //events_.resize( unique(variantList_.begin(),variantList_.end()) - variantList_.begin() );  // remove duplicates
    for (eagle::genome::EventIterator it = variantList_.begin(); it != variantList_.end(); ++it)
    {   // compiler should optimize this away in Release mode
        EAGLE_DEBUG( 0, "... " << it->getStructuralVariant() );
    }
    clog << "Sorted " << variantList_.size() << " event" << (variantList_.size() - 1 ? "s " : " ")
              << "in " << eagle::common::displayTime(clock() - start, timeProcessing) << endl;

    variantList_.chromosomeNameCheck( genome_.allContigNames() );

    // Check if some translocations appear multiple times (for example if the user specifies the same input file multiple times).
    // Warn and remove duplicated translocations
    variantList_.removeDuplicatedTranslocations();

    start = clock();
    clog << "Pairing variant list..." << endl;
    variantList_.pairing();
    clog << "Paired " << variantList_.size() << " event" << (variantList_.size() - 1 ? "s " : " ")
              << "in " << eagle::common::displayTime(clock() - start, timeProcessing) << endl;

    vector<string> shortContigIds;
    for (unsigned int dir=0; dir<2; ++dir)
    {
        for (unsigned int contigNum=0; contigNum<forwardStartingPoints.size(); ++contigNum)
        {
            const string&       contigId     = contigIds    [contigNum];
            const unsigned long contigLength = contigLengths[contigNum];
            const string contigNameRemainder = "";//contigName.substr( contigId.size() );
            eagle::model::Direction direction = dir==0 ? eagle::model::Direction::FWD : eagle::model::Direction::REV;
            unsigned int n = direction.isFwd() ? forwardStartingPoints[contigNum] : reverseStartingPoints[contigNum];

            EAGLE_DEBUG(5, "[starting points] Started processing " << n << " allele(s) from " << ((direction.isFwd())?"beginning":"end") << " of chromosome " << contigId);

            sample_.clear();
            sample_.resize(n);
//        for(unsigned int i=0; i<n; i++)
            for(eagle::genome::iterator contigOut=sample_.begin(); contigOut != sample_.end(); ++contigOut)
            {
                unsigned int i = static_cast<unsigned int>(contigOut - sample_.begin());
                contigOut->name( (boost::format("%s%s_allele%d") % prefixToAdd_ % contigId % (i+1) % ((direction.isFwd()) ? "" : "_rev") ).str()
                                 , contigNameRemainder );

                eagle::model::Locus finalPosition(contigId,(direction.isFwd())?contigLength:0);
                start = clock();
                unsigned int eventCount;
                clog << "Applying structural variants to '" << contigId << "' allele " << (i+1) << "..." << endl;
                //                EventIterator startPosition( (direction.isRev())?(bounds.second-1):bounds.first );
                eagle::genome::EventIterator startPosition = variantList_.findFirstEventForChromosome( contigId, direction );

                vector<bool> appliedStructuralVariants;
                eventCount = process( i+1, startPosition, *contigOut, direction, finalPosition
                                      , appliedStructuralVariants
                                    );

                // Save bit vector of applied structural variants
                if (appliedStructuralVariants.size() > 1 || contigLength > 1000000)
                {
/*
                    ofstream os( (outputReference_ / (contigOut->id()+".struct")).c_str() );
                    os.write( reinterpret_cast<char*>(&(appliedStructuralVariants[0])), (appliedStructuralVariants.size() + 7)/8 );
*/
                    ofstream os2( (outputReference_ / (contigOut->id()+".struct")).c_str() );
                    boost::archive::binary_oarchive oar(os2);
                    oar << appliedStructuralVariants;
                }
                else
                {
                    shortContigIds.push_back( contigOut->id() );
                }

                clog << "Applied " << eventCount << " event" << (eventCount==1 ? " " : "s ")
                     << "to a " << variantList_.ploidy().label(n) << " chromosome "
                     << "in " << eagle::common::displayTime(clock() - start, timeProcessing) << endl;

                // Remove the starting point that we just processed
                --(forwardStartingPoints[contigNum]);
                // Remove the final position that we just reached from lists of startingPoints
                vector<unsigned int >& vectorToSearch = (0 == finalPosition.pos())?forwardStartingPoints:reverseStartingPoints;
                for (unsigned int j=0; j<vectorToSearch.size(); ++j)
                {
                    if (contigIds[j] == finalPosition.chr()) {
                        EAGLE_DEBUG(5, "[starting points] Deleting starting point at " << ((0 == finalPosition.pos())?"beginning":"end") << " of chromosome " << finalPosition.chr());
                        if ( vectorToSearch[j] <= 0 ) {
                            EAGLE_ERROR("[starting points] Trying to delete already-deleted starting points at " + string((1 == finalPosition.pos())?"beginning":"end") + " of chromosome " + finalPosition.chr() + "\n You probably have a translocation landing inside a deletion\n  Vote for ticket SAGE-174 if you want this automatically detected :p" );
                        }
                        --(vectorToSearch[j]);
                        break;
                    }
                }

                start = clock();
            }
/*
            clog << "Saving " << sample_.fileCount() << " sample genome" << (sample_.fileCount() - 1 ? "s " : " ") << "..." << endl;
            sample_.save();
            clog << "Saved " <<  sample_.contigCount() << " contig" << (sample_.contigCount() - 1 ? "s " : " ")
                 << "in " << eagle::common::displayTime(clock() - start, timeIO) << endl;
*/
        }
    }

    if (shortContigIds.size())
    {
        ofstream os( (outputReference_ / "shortContigs").c_str(), ios_base::out | ios_base::app );
        BOOST_FOREACH( const string& contigId, shortContigIds )
        {
            os << contigId << endl;
        }
    }
    start = clock();
//    sample_.saveMetadata();
    clog << "Saved " GENOMESIZE_XML
              << " in " << eagle::common::displayTime(clock() - start, timeIO) << endl;

    start = clock();
    unsigned int totalVariantCount(0);
    // for (unsigned int i=0; i<variantList_.ploidy().max(); i++)
    for (unsigned int i=0; i<1; i++)
    {
        clog << "Saving " << variantList_.outputFile(i) << "..." << endl;
//        unsigned int variantCount = variantList_.save(i);
        //clog << "..saved " << variantCount << " event" << (variantCount - 1 ? "s " : " ") << endl;
//        totalVariantCount += variantCount;
    }
    clog << "Saved canonical vcf (" << totalVariantCount << " event" << (totalVariantCount - 1 ? "s" : "")
              << ") in " << eagle::common::displayTime(clock() - start, timeIO) << endl;

    clog << "+ Total processing time: " << eagle::common::displayTime(timeProcessing) << endl;
    clog << "+ Total I/O time: " << eagle::common::displayTime(timeIO) << endl;
    variantList_.check();
}

#endif //ifndef DISTRIBUTED_GENOME_MUTATOR



unsigned int GenomeMutator::process(unsigned int num,
                                    eagle::genome::EventIterator startPosition,
                                    eagle::model::Contig& contigOut,
                                    eagle::model::Direction direction,
                                    eagle::model::Locus& finalPosition
#ifdef DISTRIBUTED_GENOME_MUTATOR
                                    , vector<bool>& appliedStructuralVariants
#endif //ifdef DISTRIBUTED_GENOME_MUTATOR
                                    )
{
    eagle::genome::EventIterator lastEvent( startPosition);
    eagle::genome::EventIterator event( startPosition );
    unsigned int eventCount = 0;
//    vector< eagle::genome::EventIterator > processedDeletions;

    // Update metadata DEST field with the starting position on the allele
    event->metadata_.addInfoValue( "DEST", (boost::format("%s:0") % contigOut.id()).str() );

    while(true)
    {
        // The current event is assumed to have been applied => Skip all events that happen at the same position
        bool skip;
        do {
            event += direction.offset();
            EAGLE_DEBUG(5, "Trying... " << *event);

            skip = false;
            // Check for variants with GT field
            if( event->metadata_.hasData() && event->metadata_.getData("GT").size() > 0)
            {
                unsigned int altGtIndex = event->getVariant().altGtIndex_;
                eagle::model::Genotype GT( 1, altGtIndex );
                istringstream strm( event->metadata_.getData("GT")[0] );
                strm >> GT;
                if (GT.find( static_cast<int>(num) ) != GT.end())
                {
                    // Don't skip the ones demanded to be processed by the input VCF file

                    // But skip already processed variants
                    if (event->allele_.find(num) != event->allele_.end()) { skip = true; }
                }
                else
                {
                    EAGLE_DEBUG(5, "Skipping... " << *event);
                    skip = true;
                }
            }
            else
            {
                // Skip already processed variants
                if (!event->allele_.isHomozygousRef()) { skip = true; }
            }
            // Skip variants defined for the opposite direction
            if (!direction.sameAs( event->getVariant().adjacency.first.dir )) { skip = true; }
            // Don't skip begin/end markers
            if (!event->incoming().defined()) { skip = false; }
            // Skip variants starting at same location
            if (event->getVariant().adjacency.first.hasSameLocus( lastEvent->getVariant().adjacency.second) ) { skip = true; }
/*
            // Skip variants that would land in a deleted segment
            if (!skip)
            {
                model::Breakend dest = event->getVariant().adjacency.second;
                BOOST_FOREACH( const eagle::genome::EventIterator& del, processedDeletions )
                {
                    model::Breakend delBegin = del->getVariant().adjacency.first;
                    model::Breakend delEnd   = del->getVariant().adjacency.second;
                    if (delBegin < dest && dest < delEnd)
                    {
                        skip = true;
                        break;
                    }
                }
            }
*/

            if (skip)
            {
                EAGLE_DEBUG(5, "Skipping... " << *event);
            }
#ifdef DISTRIBUTED_GENOME_MUTATOR
            appliedStructuralVariants.push_back( !skip );
#endif //ifdef DISTRIBUTED_GENOME_MUTATOR
        } while (skip);

/*
        if (event->hasDeletion() && !event->hasTranslocation())
        {
            processedDeletions.push_back( event );
        }
*/

        // Apply variant
#ifdef DISTRIBUTED_GENOME_MUTATOR
        event->apply2(contigOut, lastEvent, make_pair(genome_.begin(),genome_.end()), direction);
#else
        event->apply(contigOut, lastEvent, make_pair(genome_.begin(),genome_.end()), direction);
#endif //ifdef DISTRIBUTED_GENOME_MUTATOR
        ++eventCount;

        // Mark as processed
        // Don't mark begin/end markers
        if( !event->allele_.set(num) )
        {
            EAGLE_WARNING("Tried to overwrite the following event");
            EAGLE_WARNING_CONT("      " << *event);
//            EAGLE_ERROR("The last warning is known to happen when a translocation break-end jumps inside a deletion. You need to fix this in the VCF file.");
        }

        // Mark paired event as processed
        int pairedEvent = event->pairedEvent_;
        if (pairedEvent)
        {
            if( !(variantList_[pairedEvent].allele_.set(num)) && !(variantList_[pairedEvent].incoming().isBiDir()))
            {
                EAGLE_WARNING("Tried to overwrite the following event");
                EAGLE_WARNING_CONT("      " << variantList_[pairedEvent]);
//                EAGLE_ERROR("The last warning is known to happen when a translocation break-end jumps inside a deletion. You need to fix this in the VCF file.");
            }
        } else {
            break;    // If we were processing a begin/end marker, leave the loop
        }

        lastEvent = event;
        // Only set a new direction if the event wasn't bi-directional.
        if (!event->outgoing().isBiDir())
        {
            direction = event->outgoing();
        }

        if (event->pairedEvent_ == (event - variantList_.begin()))
        {
            EAGLE_DEBUG(5, "Continuing to next event");
        }
        else
        {
            // Convert pointer to iterator
            event = variantList_.begin() + event->pairedEvent_;
        }
    }

    finalPosition = event->getVariant().adjacency.second;
    return eventCount;
}


} // namespace main
} // namespace eagle
