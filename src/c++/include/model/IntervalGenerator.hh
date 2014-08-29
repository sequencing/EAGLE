/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Generation of random/uniform intervals with precise boundaries
 ** Based on:
 **       [1] Jon Louis Bentley and James B. Saxe. 1980. Generating Sorted Lists of Random Numbers.
 **           ACM Trans. Math. Softw. 6, 3 (September 1980), 359-364.
 **           DOI=10.1145/355900.355907 http://doi.acm.org/10.1145/355900.355907
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MODEL_INTERVAL_GENERATOR_HH
#define EAGLE_MODEL_INTERVAL_GENERATOR_HH

#include <vector>
#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "model/FragmentLengthDist.hh"


namespace eagle
{
namespace model
{

class RandomSequenceGenerator
{
public:
    RandomSequenceGenerator(double sampleCount, double intervalSize) : remainingSampleCount_(sampleCount), intervalSize_(intervalSize), curMax_(1.0) {}
    bool hasFinished();
    unsigned long getNext();
    double getNextAsDouble();

private:
    double remainingSampleCount_;
    double intervalSize_;
    double curMax_;
};

class IntervalGenerator
{
public:
    IntervalGenerator( const std::vector<unsigned long>& contigLengths, const unsigned long readCount, const unsigned int minFragmentLength, const unsigned int medianFragmentLength, const unsigned int maxFragmentLength, const bool verbose)
        : contigLengths_        ( contigLengths )
        , minFragmentLength_    ( minFragmentLength )
        , medianFragmentLength_ ( medianFragmentLength )
        , maxFragmentLength_    ( maxFragmentLength )
        , verbose_              ( verbose )
    {}
    virtual ~IntervalGenerator() {}
    virtual std::pair< unsigned long, unsigned int > getNext( const signed long testValue = -1 ) = 0;

protected:
    const std::vector<unsigned long> contigLengths_;
    const unsigned int minFragmentLength_;
    const unsigned int medianFragmentLength_;
    const unsigned int maxFragmentLength_;
    bool verbose_;
};

class RandomIntervalGenerator : public IntervalGenerator
{
public:
    RandomIntervalGenerator( const std::vector<unsigned long>& contigLengths, const unsigned long readCount, const unsigned int minFragmentLength, const unsigned int medianFragmentLength, const unsigned int maxFragmentLength, const bool verbose=true)
        : IntervalGenerator( contigLengths, readCount, minFragmentLength, medianFragmentLength, maxFragmentLength, verbose)
        , intervalsPerNormalPosition_( maxFragmentLength_ - minFragmentLength_ + 1)
        , contigIntervalInfo_        ( )
        , randomSeq_                 ( readCount, getIntervalCount() )
        , currentContigIntervalInfo_ ( contigIntervalInfo_.begin() )
    {
        assert( maxFragmentLength_ >= minFragmentLength_ );
    }

    virtual ~RandomIntervalGenerator() {}
    virtual std::pair< unsigned long, unsigned int > getNext( const signed long testValue = -1 );

private:
    unsigned long getIntervalCount();

    struct ContigIntervalInfo {
        unsigned long firstInterval;
        unsigned long lastInterval;
        unsigned long firstSmallerInterval;
        unsigned long firstIntervalGlobalPos;
        unsigned long firstSmallerIntervalGlobalPos;
    };

    const unsigned int intervalsPerNormalPosition_;
    std::vector<ContigIntervalInfo> contigIntervalInfo_;
    RandomSequenceGenerator randomSeq_;
    std::vector<ContigIntervalInfo>::iterator currentContigIntervalInfo_;
};

class RandomIntervalGeneratorUsingIntervalLengthDistribution : public IntervalGenerator
{
public:
    RandomIntervalGeneratorUsingIntervalLengthDistribution( const std::vector<unsigned long>& contigLengths, const unsigned long readCount, const boost::filesystem::path& templateLengthTableFile, const bool verbose=true)
        : IntervalGenerator( contigLengths, readCount, 0, 0, 0, verbose)
        , currentPos_( 0 )
        , accumulatedProbasUntilCurrentPos_( 0 )
        , fragmentLengthDist_( templateLengthTableFile )
        , randomSeq_( readCount, getTotalIntervalsProba() )
    {
    }

    virtual ~RandomIntervalGeneratorUsingIntervalLengthDistribution() {}
    virtual std::pair< unsigned long, unsigned int > getNext( const signed long testValue = -1 );

private:
    double getIntervalsProbaAtPos( unsigned long globalPos );
    double getTotalIntervalsProba();

    std::vector<float> probas_;
    unsigned long currentPos_;
    double accumulatedProbasUntilCurrentPos_;
    FragmentLengthDist fragmentLengthDist_;
    RandomSequenceGenerator randomSeq_;
};

class UniformIntervalGenerator : public IntervalGenerator
{
public:
    UniformIntervalGenerator( const std::vector<unsigned long>& contigLengths, const unsigned int medianFragmentLength, const double step, unsigned long& readCount, const bool verbose=true);
    virtual ~UniformIntervalGenerator() {}
    virtual std::pair< unsigned long, unsigned int > getNext( const signed long testValue = -1 );

private:
    double step_;
    struct ContigIntervalInfo
    {
        unsigned long firstGlobalPos;
        unsigned long lastGlobalPos;
    };
    std::vector<ContigIntervalInfo> contigIntervalInfo_;
    std::vector<ContigIntervalInfo>::iterator currentContigIntervalInfo_;
    double currentGlobalPos_;
};

class RandomIntervalGeneratorFromProbabilityMatrix : public IntervalGenerator
{
public:
    RandomIntervalGeneratorFromProbabilityMatrix( const FragmentLengthProbabilityMatrix& fragmentLengthProbabilityMatrix, const std::vector<unsigned long>& contigLengths, unsigned long& readCount, const bool verbose=true)
        : IntervalGenerator( contigLengths, readCount, 0, 0, 0, false )
        , fragmentLengthProbabilityMatrix_( fragmentLengthProbabilityMatrix )
        , fragmentLengthDist_min( fragmentLengthProbabilityMatrix.fragmentLengthDist_.min() )
        , fragmentLengthDist_size( fragmentLengthProbabilityMatrix.fragmentLengthDist_.size() )
        , intervalSize_( fragmentLengthProbabilityMatrix.sum_P_FL_pos )
        , curMax_( 1.0 )
        , remainingSampleCount_( readCount )
        , lastChoice( 0.0 )
        , lastChoiceIndex( 0 )
        , lastChoiceRemainder( 0.0 )
        {}
    virtual ~RandomIntervalGeneratorFromProbabilityMatrix() {}
    virtual std::pair< unsigned long, unsigned int > getNext( const signed long testValue = -1 );

private:
    model::FragmentLengthProbabilityMatrix fragmentLengthProbabilityMatrix_;

    const unsigned int fragmentLengthDist_min;
    const unsigned int fragmentLengthDist_size;
    const double intervalSize_;
    double curMax_;
    unsigned int remainingSampleCount_;
    double lastChoice;
    unsigned int lastChoiceIndex;
    double lastChoiceRemainder;
};



} // namespace model
} // namespace eagle

#endif //EAGLE_MODEL_INTERVAL_GENERATOR_HH
