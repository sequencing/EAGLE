/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Top level component to induce variants in a reference.
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MAIN_VCF_COMPARATOR_HH
#define EAGLE_MAIN_VCF_COMPARATOR_HH

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "genome/VariantList.hh"
#include "main/VcfComparatorOptions.hh"

namespace bfs = boost::filesystem;


namespace eagle
{
namespace main
{

class VcfComparator
{
public:
    VcfComparator( const VcfComparatorOptions& options );
    void run();

private:
    const VcfComparatorOptions &options_;
    eagle::genome::VariantList simulatedVariantList_;
    eagle::genome::VariantList calledVariantList_;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_VCF_COMPARATOR_HH
