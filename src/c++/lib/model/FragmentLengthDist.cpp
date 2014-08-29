/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include <numeric>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "model/FragmentLengthDist.hh"

using namespace std;


namespace eagle
{
namespace model
{


FragmentLengthDist::FragmentLengthDist( const boost::filesystem::path& filename )
    : min_(0)
    , max_(0)
{
    // Parse template length table
    io::DsvReader tsvReader( filename );
    vector<string> tokens;
    while ( tsvReader.getNextLineFields<'\t'>(tokens) )
    {
        if ( tokens.size() == 0 ) { continue; }
        assert ( tokens.size() == 2 && "There should be 2 entries per line" );
        vector<unsigned long> values;
        unsigned int templateLength;
        double templateLengthCount;
        try
        {
            templateLength = boost::lexical_cast<unsigned long,std::string>( tokens[0] );
            templateLengthCount = boost::lexical_cast<double,std::string>  ( tokens[1] );
        }
        catch (const boost::bad_lexical_cast &e)
        {
            EAGLE_ERROR("Error while reading template length table: a numerical field seems to contain non-numerical characters");
        }

        if (min_ == 0)
        {
            min_ = templateLength;
        }
        if (templateLength < min_)
        {
            EAGLE_ERROR("Error: the template length table should be sorted by template length (first column)");
        }
        if (templateLength > max_)
        {
            max_ = templateLength;
            templateLengthDist_.resize( max_ - min_ + 1 );
        }
        assert( templateLength >= min_ );
        assert( templateLength - min_ < templateLengthDist_.size() );
        templateLengthDist_[ templateLength - min_ ] = templateLengthCount;
    }

    // Normalise table
    double sum = std::accumulate( templateLengthDist_.begin(), templateLengthDist_.end(), 0.0 );
    if (sum == 0)
    {
        EAGLE_ERROR("Error: the template length table counts (column 2) adds up to zero");
    }
    BOOST_FOREACH( double& val, templateLengthDist_ )
    {
        val /= sum;
    }
}


void FragmentLengthProbabilityMatrix::phase1( const unsigned long chrLength, const FragmentLengthDist& fragmentLengthDist )
{
    // For each position, use GC and CNV to establish the probability of a fragment covering this pos:
    // P(pos) = GC_trans[GC%(pos)] * CNV(pos)
    // Both GC% and CNV are faster to calculate if we keep the value for the previous position
    GcPercentage gcPercentage;
    CNV cnv;
    GcTrans gcTrans;
    vector<float> P( chrLength );
    for (unsigned int globalPos = 0; globalPos < chrLength; ++globalPos)
    {
        double gcPercentage_pos = gcPercentage[ globalPos ];
        double gcMultiplier_pos = gcTrans[ gcPercentage_pos ];
        double cnvMultiplier_pos = cnv[ globalPos ];
        double P_pos = gcMultiplier_pos * cnvMultiplier_pos;
        P[globalPos] = P_pos;
    }

#define COUNT_ONLY
#ifdef COUNT_ONLY
    P_FL_pos->resize( fragmentLengthDist.size() );
#else
    P_FL_pos->resize( fragmentLengthDist.size() * chrLength );
#endif  //ifdef COUNT_ONLY
    sum_P_FL_pos = 0.0;
    for (unsigned int globalPos = 0; globalPos < chrLength; ++globalPos)
    {
        double P_FL_pos_tmp = 1.0;
        unsigned int P_FL_pos_index = 0;//globalPos * fragmentLengthDist.size();
        unsigned int FLmax = min<unsigned int>( fragmentLengthDist.max(), chrLength - globalPos );
        if (globalPos%100000 == 0)
        {
            cout << "calculating for globalPos=" << globalPos << ", FLmax=" << FLmax << endl;
        }
        for (unsigned int FL = fragmentLengthDist.min(); FL <= FLmax; ++FL, ++P_FL_pos_index)
        {
/*
            P_FL_pos[FL][globalPos] = 1.0;
            for (unsigned int x = 0; x < FL; ++x)
            {
                P_FL_pos[FL][globalPos] *= P[globalPos + x];
            }
*/
            P_FL_pos_tmp *= P[globalPos + FL];
            P_FL_pos->operator[](P_FL_pos_index) = P_FL_pos_tmp;
        }
        P_FL_pos_index = 0;//globalPos * fragmentLengthDist.size();
        for (unsigned int FL = fragmentLengthDist.min(); FL <= FLmax; ++FL, ++P_FL_pos_index)
        {
            double tmp = P_FL_pos->operator[](P_FL_pos_index) * fragmentLengthDist[FL];
            sum_P_FL_pos += tmp;
        }
//        cout << "sum_P_FL_pos=" << sum_P_FL_pos << endl;
    }
    cout << "Final enumeration: " << sum_P_FL_pos << endl;
    exit(0);
/*
    // Validity check
    {
        double sum = 0;
        for (unsigned int i=0; i<P_FL_pos->size(); ++i)
        {
            sum += P_FL_pos[i];
        }
        cout << "sum=" << sum << endl;
    }
*/
}


} // namespace model
} // namespace eagle
