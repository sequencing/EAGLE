/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MODEL_FRAGMENT_LENGTH_DIST_HH
#define EAGLE_MODEL_FRAGMENT_LENGTH_DIST_HH

#include <vector>
#include <boost/filesystem.hpp>

#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "io/Text.hh"


namespace eagle
{
namespace model
{


class GcPercentage
{
public:
    double operator[]( const unsigned int globalPos ) const
    {
        return 0.5;
    }
};

class CNV
{
public:
    double operator[]( const unsigned int globalPos ) const
    {
        return 1.0;
    }
};

class GcTrans
{
public:
    double operator[]( const double gcPercentage ) const
    {
        return 1.0;
    }
};

class BlockMatrix
{
public:
    double  operator[]( const unsigned long pos ) const { return val0_; }
    double& operator[]( const unsigned long pos )       { return val0_; }
    unsigned long size() const { return 0; }

private:
    double val0_;
};

class FragmentLengthDist
{
public:
    FragmentLengthDist( const boost::filesystem::path& filename );
    unsigned int min() const { return min_; }
    unsigned int max() const { return max_; }
    unsigned int size() const { return max() - min() + 1; }
    double operator[]( const unsigned int templateLength ) const
    {
        assert( templateLength >= min_ );
        assert( templateLength - min_ < templateLengthDist_.size() );
        return templateLengthDist_[ templateLength - min_ ];
    }

private:
    unsigned int min_;
    unsigned int max_;
    std::vector<double> templateLengthDist_;
};

class FragmentLengthProbabilityMatrix
{
public:
    FragmentLengthProbabilityMatrix( const unsigned long chrLength, const boost::filesystem::path& templateLengthTableFilename )
        : fragmentLengthDist_( templateLengthTableFilename )
        , P_FL_pos( new std::vector< double >() )
    {
        phase1( chrLength, fragmentLengthDist_ );
//        assert( P_FL_pos.get() != 0 );
    }
/*
    FragmentLengthProbabilityMatrix( const boost::filesystem::path& blockCoverageFilename, const boost::filesystem::path& templateLengthTableFilename )
        : fragmentLengthDist_( templateLengthTableFilename )
        , blockCoverage_( blockCoverageFilename )
    {
//        phase1( chrLength, fragmentLengthDist_ );
//        assert( P_FL_pos.get() != 0 );
    }
*/
    FragmentLengthProbabilityMatrix( const FragmentLengthProbabilityMatrix& obj )
        : fragmentLengthDist_( obj.fragmentLengthDist_ )
        , sum_P_FL_pos( obj.sum_P_FL_pos )
        , P_FL_pos( obj.P_FL_pos )
    {
//        assert( P_FL_pos.get() != 0 );
    }

    const std::vector< double >& getProbabilities() { return *P_FL_pos; };
//    const BlockMatrix& getProbabilities() { return P_FL_pos; };

    FragmentLengthDist fragmentLengthDist_;
//    BlockCoverage blockCoverage_;
    double sum_P_FL_pos;
    boost::shared_ptr< std::vector< double > > P_FL_pos;
//    BlockMatrix P_FL_pos;

private:
    void phase1( const unsigned long chrLength, const FragmentLengthDist& fragmentLengthDist );
};


} // namespace model
} // namespace eagle

#endif //EAGLE_MODEL_FRAGMENT_LENGTH_DIST_HH
