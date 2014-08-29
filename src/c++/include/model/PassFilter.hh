/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Class used to determine whether a cluster will Pass Filter (PF) or not
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MODEL_PASS_FILTER_HH
#define EAGLE_MODEL_PASS_FILTER_HH

#include <string>


namespace eagle
{
namespace model
{

class PassFilter
{
public:
    static bool isBclClusterPassingFilter( const char *bclCluster, const unsigned int clusterLength )
    {
        unsigned int Ncount = 0;
        for (unsigned int i=0; i<clusterLength; ++i)
        {
            Ncount += (bclCluster[i] == 0)?1:0;
        }
        return (Ncount < 64);
    }

    static bool isSequencePassingFilter( const std::string& seq )
    {
        unsigned int Ncount = std::count( seq.begin(), seq.end(), 'N' );
        return (Ncount < 64);
    }
};


} // namespace model
} // namespace eagle

#endif // EAGLE_MODEL_PASS_FILTER_HH
