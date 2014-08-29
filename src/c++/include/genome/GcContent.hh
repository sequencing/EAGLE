/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Set of classes to calculate and use GC content
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_GENOME_GC_CONTENT_HH
#define EAGLE_GENOME_GC_CONTENT_HH

#include <utility>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "model/Fragment.hh"
#include "model/Nucleotides.hh"


using namespace std;


namespace eagle
{
namespace genome
{


class GcCoverageFit
{
public:
    GcCoverageFit( const boost::filesystem::path& gcCoverageFitFilename, const boost::filesystem::path& sampleGenomeDir/*, boost::shared_ptr<boost::mt19937> randomGen*/ );
    double averageMultiplier();
    bool needsDiscarding( const model::Fragment& fragment );
    bool needsDiscarding( const double gcContent );

private:
    bool isActive_;
    std::vector<double> gcContentValues_;
    std::vector<double> coverageMultiplierValues_;
//    boost::shared_ptr<boost::mt19937>& randomGen_;
    eagle::model::IUPAC baseConverter_;
    double averageMultiplier_;

    void parseGcCoverageFitFile( const boost::filesystem::path& filename );
    double getInterpolatedCoverageMultiplierForGcContent( const double gcContent );
};



} // namespace genome
} // namespace eagle

#endif // EAGLE_GENOME_GC_CONTENT_HH
