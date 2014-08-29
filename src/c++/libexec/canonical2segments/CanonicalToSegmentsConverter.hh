/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MAIN_CANONICAL_TO_SEGMENTS_CONVERTER_HH
#define EAGLE_MAIN_CANONICAL_TO_SEGMENTS_CONVERTER_HH

#include "io/RunInfo.hh"
#include "genome/EnrichedFragment.hh"
#include "genome/ReadCluster.hh"
#include "CanonicalToSegmentsConverterOptions.hh"


namespace eagle
{
namespace main
{


class CanonicalToSegmentsConverter
{
public:
    CanonicalToSegmentsConverter (const CanonicalToSegmentsConverterOptions &options);
    void run();

private:
    const CanonicalToSegmentsConverterOptions &options_;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_CANONICAL_TO_SEGMENTS_CONVERTER_HH
