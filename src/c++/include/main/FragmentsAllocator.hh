/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_FRAGMENTS_ALLOCATOR_HH
#define EAGLE_FRAGMENTS_ALLOCATOR_HH

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "model/Fragment.hh"
#include "model/IntervalGenerator.hh"
#include "genome/GcContent.hh"
#include "FragmentsAllocatorOptions.hh"


namespace eagle
{
namespace main
{

class FragmentsAllocator
{
public:
    FragmentsAllocator( const FragmentsAllocatorOptions &options );
    void run();

private:
    void setRandomSeed();
    model::FragmentWithAllocationMetadata getNextFragment( boost::shared_ptr<eagle::model::IntervalGenerator>& randomInterval,
                                                           boost::shared_ptr<eagle::model::MultiFragmentFilesReader>& multiFragmentFilesReader,
                                                           const unsigned long fragmentNum,
                                                           const unsigned long fragmentCount );

    const FragmentsAllocatorOptions &options_;
//    boost::shared_ptr< boost::mt19937 > randomGen_;
    genome::GcCoverageFit gcCoverageFit_;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_FRAGMENTS_ALLOCATOR_HH
