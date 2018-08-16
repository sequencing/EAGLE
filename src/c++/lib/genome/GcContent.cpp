/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include <boost/filesystem.hpp>
#include "common/Logger.hh"
#include "genome/SharedFastaReference.hh"
#include "io/Text.hh"

#include "genome/GcContent.hh"


using namespace std;


namespace eagle
{
namespace genome
{


GcCoverageFit::GcCoverageFit( const boost::filesystem::path& gcCoverageFitFilename, const boost::filesystem::path& sampleGenomeDir/*, boost::shared_ptr<boost::mt19937> randomGen*/ )
    : isActive_( !gcCoverageFitFilename.empty() )
                            //    , randomGen_( randomGen )
{
    if (isActive_)
    {
        parseGcCoverageFitFile( gcCoverageFitFilename );
        // already done somewhere else: genome::SharedFastaReference::init( sampleGenomeDir );
    }
}

void GcCoverageFit::parseGcCoverageFitFile( const boost::filesystem::path& filename )
{
    io::DsvReader tsvReader( filename );
    vector<string> tokens;
    double maxCovMult = 0;
    double weightedSum = 0;
    double totalWeight = 0;

    while( tsvReader.getNextLineFields<'\t'>(tokens) )
    {
        assert( (tokens.size() == 2 || tokens.size() == 3) && "GC_coverage_fit table should have 2 or 3 entries per line: gc%, normalised coverage (coverage multiplier) and an optional window count in ref" );
        vector<double> values;
        try
        {
            std::transform( tokens.begin(), tokens.end(), std::back_inserter(values), boost::bind( &boost::lexical_cast<double,std::string>, _1) );
        }
        catch (const boost::bad_lexical_cast &e)
        {
            EAGLE_ERROR("Error while reading GC_coverage_fit table: a field seems to contain non-numerical characters");
        }
        assert( values.size() == tokens.size() );
        assert( values[0] >= 0 && "Error while reading GC_coverage_fit table: Column 1 values should be >= 0" );
        assert( values[0] <= 100 && "Error while reading GC_coverage_fit table: Column 1 values should be <= 100" );

        gcContentValues_.push_back( values[0]/100 );
        coverageMultiplierValues_.push_back( values[1] );
        if (values.size() >= 3)
        {
            weightedSum += values[1] * values[2];
            totalWeight += values[2];
        }
        maxCovMult = max<double>( maxCovMult, values[1] );
    }

    assert( gcContentValues_.size() == coverageMultiplierValues_.size() );
    assert( gcContentValues_.size() > 0  && "No valid entries found in GC_coverage_fit table file");

    // Normalise
    BOOST_FOREACH( double& covMult, coverageMultiplierValues_ )
    {
        covMult /= maxCovMult;
    }
    weightedSum /= maxCovMult;

    averageMultiplier_ = totalWeight ? ( weightedSum / totalWeight ) : getInterpolatedCoverageMultiplierForGcContent( 0.44 );
}

double GcCoverageFit::averageMultiplier()
{
    if (isActive_)
    {
        return averageMultiplier_;
    }
    else
    {
        return 1.0;
    }
}

bool GcCoverageFit::needsDiscarding( const model::Fragment& fragment )
{
    if (!isActive_)
    {
        return false;
    }

//#define EXPERIMENT_SPECIAL_GC_BIAS_AT_CYCLE_4
#ifdef EXPERIMENT_SPECIAL_GC_BIAS_AT_CYCLE_4
    if (fragment.fragmentLength_ > 4)
    {
        unsigned long offset = 3;
        bool overlapContigBoundary;
        char base = genome::SharedFastaReference::get()->get( fragment.startPos_, offset, overlapContigBoundary );
        base = baseConverter_.norm( base );
        base = toupper( base );
        if (base=='A' || base=='T')
        {
            double randomValue = (double)random() / (double)RAND_MAX;
            if (randomValue < .25)
                return true; //discard
        }

        offset = 4;
        base = genome::SharedFastaReference::get()->get( fragment.startPos_, offset, overlapContigBoundary );
        base = baseConverter_.norm( base );
        base = toupper( base );
        if (base=='A' || base=='T')
        {
            double randomValue = (double)random() / (double)RAND_MAX;
            if (randomValue < .1)
                return true; //discard
        }
    }
#endif

    unsigned int gcCount=0, acgtCount=0;
    for (unsigned long offset = 0; offset < fragment.fragmentLength_; ++offset)
    {
        // Only consider the first and last 150 bases, to make the GC% values less average and to match a bit more closely Firebrand's GC plots, which are calculated using 150bp windows
        if (offset >= 150 && offset < fragment.fragmentLength_-150)
        {
            continue;
        }

        bool overlapContigBoundary;
        char base = genome::SharedFastaReference::get()->get( fragment.startPos_, offset, overlapContigBoundary );
        if (overlapContigBoundary)
        {
            // safety check
            return true;
        }
        base = baseConverter_.norm( base );
        base = toupper( base );

        switch (base)
        {
        case 'C':
        case 'G':
            ++gcCount;
        case 'A':
        case 'T':
            ++acgtCount;
            break;
        default:
            return needsDiscarding( averageMultiplier() ); // if, for example, any 'N' appears, we keep the fragment
        }
    }

    // special case
    if (acgtCount == 0)
    {
        return needsDiscarding( averageMultiplier() );
    }

    double gcContent = ((double)gcCount) / acgtCount;
    return needsDiscarding( gcContent );
}

bool GcCoverageFit::needsDiscarding( const double gcContent )
{
    double covMult = getInterpolatedCoverageMultiplierForGcContent( gcContent );

    if (covMult == 1.0)
    {
        return false;
    }

//    double randomValue = (double)randomGen_->operator()() / (double)randomGen_->max();
    double randomValue = (double)random() / (double)RAND_MAX;
    return (randomValue > covMult);
}

double GcCoverageFit::getInterpolatedCoverageMultiplierForGcContent( const double gcContent )
{
    std::vector<double>::iterator itr = std::lower_bound( gcContentValues_.begin(), gcContentValues_.end(), gcContent );
    if (itr == gcContentValues_.begin())
    {
        return coverageMultiplierValues_.front();
    }
    if (itr == gcContentValues_.end())
    {
        return coverageMultiplierValues_.back();
    }
    unsigned int index = itr - gcContentValues_.begin();
    if (*itr == gcContent)
    {
        return coverageMultiplierValues_[index];
    }
    // Last case: we need interpolation
    double interpolationFactor = (gcContent - gcContentValues_[index-1]) / (gcContentValues_[index] - gcContentValues_[index-1]);
    double interpolatedResult = coverageMultiplierValues_[index-1] + interpolationFactor * (coverageMultiplierValues_[index] - coverageMultiplierValues_[index-1]);
    return interpolatedResult;
}


} // namespace genome
} // namespace eagle
