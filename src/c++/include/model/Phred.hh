/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Pre-computed Phred quality scores
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MODEL_PHRED_HH
#define EAGLE_MODEL_PHRED_HH

#include <math.h>
#include <algorithm>
#include <vector>
#include <cassert>
#include "common/Exceptions.hh"


namespace eagle
{
namespace model
{

class Phred
{
public:
    static const unsigned int QUALITY_MAX = 50;

//    Phred() {}
    static inline void initCheck()
    {
        if (qualityToProbability_.size() == 0)
        {
            qualityToProbability_.resize(QUALITY_MAX+1);
            for (unsigned int q=0; q<=QUALITY_MAX; ++q)
            {
                qualityToProbability_[q] = pow(10, -(double(q) / 10));
            }
        }
    }

    static double qualToProb( const unsigned int qual )
    {
        initCheck();
        if (qual > QUALITY_MAX)
        {
            BOOST_THROW_EXCEPTION( eagle::common::EagleException( 0, "Phred quality is higher than allowed max") );
        }
        double prob = qualityToProbability_[qual];
        return prob;
    }

    static unsigned int probToQual( const double prob )
    {
        initCheck();
        std::vector< double >::iterator it = 
            std::lower_bound( qualityToProbability_.begin(), qualityToProbability_.end(), prob, std::greater<double>() );
        unsigned int qual = it - qualityToProbability_.begin();
        return qual;
    }

private:
    static std::vector< double > qualityToProbability_;
};



} // namespace model
} // namespace eagle

#endif // EAGLE_MODEL_PHRED_HH
