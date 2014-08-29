/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include <boost/assign.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/integer.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <fstream>
#include <cmath>

#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "main/VcfComparator.hh"

using namespace std;


namespace eagle
{
namespace main
{


VcfComparator::VcfComparator( const VcfComparatorOptions &options )
    : options_              ( options )
    , simulatedVariantList_ ( options.simulatedVariants, "", 1 )
    , calledVariantList_    ( options.calledVariants, "", 1 )
{
}


void VcfComparator::run()
{
#ifdef EAGLE_DEBUG_MODE
    unsigned long timeProcessing = 0;
    unsigned long timeIO = 0;
#endif //ifdef EAGLE_DEBUG_MODE
    clock_t start;

    // Load simulated variants
    start = clock();
    EAGLE_DEBUG( 0, (boost::format("Loading %d simulated variant list(s)...") % simulatedVariantList_.fileCount()).str() );
    simulatedVariantList_.load();
    EAGLE_DEBUG( 0, (boost::format("Loaded %d event(s) in ") % simulatedVariantList_.size()).str()
                 << eagle::common::displayTime(clock() - start, timeIO) );

    // Sort them
    start = clock();
    EAGLE_DEBUG( 0, "Sorting simulated variant list..." );
    std::stable_sort( simulatedVariantList_.begin(), simulatedVariantList_.end(), genome::Event::ltComparisonIncludingAltField );
    for (eagle::genome::EventIterator it = simulatedVariantList_.begin(); it != simulatedVariantList_.end(); ++it)
    {   // compiler should optimize this away in Release mode
        EAGLE_DEBUG( 0, "... " << it->getStructuralVariant() );
    }
    EAGLE_DEBUG( 0, (boost::format("Sorted %d event(s) in ") % simulatedVariantList_.size()).str()
                 << eagle::common::displayTime(clock() - start, timeProcessing) );


    // Load called variants
    start = clock();
    EAGLE_DEBUG( 0, (boost::format("Loading %d called variant list(s)...") % calledVariantList_.fileCount()).str() );
    calledVariantList_.load();
    EAGLE_DEBUG( 0, (boost::format("Loaded %d event(s) in ") % calledVariantList_.size()).str()
                 << eagle::common::displayTime(clock() - start, timeIO) );

    // Sort them
    start = clock();
    EAGLE_DEBUG( 0, "Sorting called variant list..." );
    std::stable_sort( calledVariantList_.begin(), calledVariantList_.end(), genome::Event::ltComparisonIncludingAltField );
    for (eagle::genome::EventIterator it = calledVariantList_.begin(); it != calledVariantList_.end(); ++it)
    {   // compiler should optimize this away in Release mode
        EAGLE_DEBUG( 0, "... " << it->getStructuralVariant() );
    }
    EAGLE_DEBUG( 0, (boost::format("Sorted %d event(s) in ") % calledVariantList_.size()).str()
                 << eagle::common::displayTime(clock() - start, timeProcessing) );


    // Go through both lists of variants and compare them
    eagle::genome::EventIterator simulatedVariantItr = simulatedVariantList_.begin();
    eagle::genome::EventIterator calledVariantItr = calledVariantList_.begin();
    while (simulatedVariantItr != simulatedVariantList_.end() || calledVariantItr != calledVariantList_.end())
    {
        int diff;
        if (simulatedVariantItr == simulatedVariantList_.end())
        {
            diff = -1;
        }
        else if (calledVariantItr == calledVariantList_.end())
        {
            diff = +1;
        }
        else if (*calledVariantItr < *simulatedVariantItr) // comparison by position
        {
            diff = -1;
        }
        else if (*simulatedVariantItr < *calledVariantItr) // comparison by position
        {
            diff = +1;
        }
        else if (*simulatedVariantItr == *calledVariantItr) // comparison by full variant info
        {
            // We need to check ahead whether all the calls at the same position are equal to the simulated variants
            unsigned int sameSimulatedVariantCount = 0;
            while ( (simulatedVariantItr+sameSimulatedVariantCount+1 != simulatedVariantList_.end()) && !(*calledVariantItr < *(simulatedVariantItr+sameSimulatedVariantCount+1)) )
            {
                ++sameSimulatedVariantCount;
            }
            unsigned int sameCalledVariantCount = 0;
            while ( (calledVariantItr+sameCalledVariantCount+1 != calledVariantList_.end()) && !(*simulatedVariantItr < *(calledVariantItr+sameCalledVariantCount+1)) )
            {
                ++sameCalledVariantCount;
            }
            if (sameSimulatedVariantCount != sameCalledVariantCount)
            {
                diff = 2;
            }
            else
            {
                bool allVariantsMatch = true;
                for (unsigned int i=0; i<sameSimulatedVariantCount; ++i)
                {
                    if (!(*(simulatedVariantItr+i) == *(calledVariantItr+i)))
                    {
                        allVariantsMatch = false;
                        break;
                    }
                }
                if (allVariantsMatch)
                {
                    diff = 0;
                }
                else
                {
                    diff = 2;
                }
            }
        }
        else // positions are here equal, but full variants differ
        {
            diff = 2;
        }

        switch (diff)
        {
        case -1:
            if (!calledVariantItr->incoming().isRev())
            {
                cout << "Variant only in called list   (FP):\t" << calledVariantItr->metadata_.strInfo() << "\t" << *calledVariantItr << "\t" << calledVariantItr->metadata_.qual << endl;
            }
            ++calledVariantItr;
            break;
        case 1:
            if (!simulatedVariantItr->incoming().isRev())
            {
                cout << "Variant only in simulated list(FN):\t" << *simulatedVariantItr << "\t" << simulatedVariantItr->metadata_.qual << endl;
            }
            ++simulatedVariantItr;
            break;
        case 0:
            cout << "Correct call(TP):\t" << calledVariantItr->metadata_.strInfo() << "\t" << *simulatedVariantItr << " == " << *calledVariantItr << "\t" << simulatedVariantItr->metadata_.qual << "\t" << calledVariantItr->metadata_.qual << "\t" << simulatedVariantItr->metadata_.filter << "\t" << calledVariantItr->metadata_.filter << endl;
            ++simulatedVariantItr;
            ++calledVariantItr;
            break;
        case 2:
            cout << "Incorrect call(FP+FN):\t" << calledVariantItr->metadata_.strInfo() << "\t" << *simulatedVariantItr << " == " << *calledVariantItr;
            while ( (simulatedVariantItr+1 != simulatedVariantList_.end()) && !(*calledVariantItr < *(simulatedVariantItr+1)) )
            {
                ++simulatedVariantItr;
                cout << " == " << *simulatedVariantItr;
            }
            while ( (calledVariantItr+1 != calledVariantList_.end()) && !(*simulatedVariantItr < *(calledVariantItr+1)) )
            {
                ++calledVariantItr;
                cout << " == " << *calledVariantItr;
            }
            ++simulatedVariantItr;
            ++calledVariantItr;
            cout << endl;
            break;
        default:
            assert(false);
        }
    }
}


} // namespace main
} // namespace eagle
