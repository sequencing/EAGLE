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

#include <numeric>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "model/IntervalGenerator.hh"

using namespace std;


namespace eagle
{
namespace model
{

// Class RandomSequenceGenerator

unsigned long RandomSequenceGenerator::getNext()
{
    return static_cast<unsigned long>( floor(getNextAsDouble()) );
}

bool RandomSequenceGenerator::hasFinished()
{
    return (remainingSampleCount_ <= 0);
}

double RandomSequenceGenerator::getNextAsDouble()
{
    if (remainingSampleCount_ <= 0)
    {
        BOOST_THROW_EXCEPTION( eagle::common::EagleException( 0, "RandomSequenceGenerator: all samples have already been generated") );
    }
    double rand0to1 = (1.0+rand())/(1.0+RAND_MAX);
    curMax_ *= exp(log(rand0to1)/remainingSampleCount_);
    remainingSampleCount_--;
    return (1-curMax_)*intervalSize_;
}


// Class RandomIntervalGenerator

pair< unsigned long, unsigned int > RandomIntervalGenerator::getNext( const signed long testValue )
{
    unsigned long intervalNum = (testValue<0)?randomSeq_.getNext():testValue;
    
    // Convert intervalNum to interval
    while ( intervalNum > currentContigIntervalInfo_->lastInterval || currentContigIntervalInfo_->lastInterval == 0xFFFFFFFFFFFFFFFF )
    {
        ++currentContigIntervalInfo_;
        if (verbose_)
        {
            clog << "Generating fragments for chromosome starting at global pos " << currentContigIntervalInfo_->firstIntervalGlobalPos << endl;
        }
    }
    assert( intervalNum >= currentContigIntervalInfo_->firstInterval );

    pair< unsigned long, unsigned int > p;
    if ( intervalNum < currentContigIntervalInfo_->firstSmallerInterval )
    {
        unsigned long localIntervalNum = intervalNum - currentContigIntervalInfo_->firstInterval;
        p.first  = currentContigIntervalInfo_->firstIntervalGlobalPos + localIntervalNum / intervalsPerNormalPosition_;
        p.second = minFragmentLength_ + localIntervalNum % intervalsPerNormalPosition_;
    }
    else
    {
//        static unsigned long lastSmallerIntervalNum = -1;
//        if ( currentContigIntervalInfo_->firstSmallerInterval != lastSmallerIntervalNum )
        {
            unsigned long localIntervalNum = intervalNum - currentContigIntervalInfo_->firstSmallerInterval;
            unsigned int intervalCountForThisPos = maxFragmentLength_ - minFragmentLength_;
            unsigned long currentPos = currentContigIntervalInfo_->firstSmallerIntervalGlobalPos;
            while (intervalCountForThisPos > 0)
            {
                if ( localIntervalNum < intervalCountForThisPos )
                {
                    p.first  = currentPos;
                    p.second = localIntervalNum + minFragmentLength_;
                    break;
                }
                localIntervalNum -= intervalCountForThisPos;
                --intervalCountForThisPos;
                ++currentPos;
            }
            assert( intervalCountForThisPos>0 && "should never happen" );
        }
    }
    return p;
}

unsigned long RandomIntervalGenerator::getIntervalCount()
{
    unsigned long intervalCount = 0;
    unsigned long globalPos = 0;
    BOOST_FOREACH( const unsigned long l, contigLengths_ )
    {
        struct ContigIntervalInfo intervalInfo;
        intervalInfo.firstInterval = intervalCount;
        intervalInfo.firstIntervalGlobalPos = globalPos;

        if ( l+1 >= maxFragmentLength_ )
        {
            unsigned long normalPositionCount = l - maxFragmentLength_ + 1;
            intervalCount += intervalsPerNormalPosition_ * normalPositionCount;
            intervalInfo.firstSmallerInterval = intervalCount;
            intervalCount += intervalsPerNormalPosition_ * (intervalsPerNormalPosition_-1) / 2;

            intervalInfo.firstSmallerIntervalGlobalPos = globalPos + normalPositionCount;
        }
        else
        {
            intervalInfo.firstSmallerIntervalGlobalPos = 0;
            EAGLE_WARNING( "Chromosome shorter than insert length" );
        }
        intervalInfo.lastInterval = intervalCount - 1;

        if (verbose_)
        {
            clog << "Chromosome at global pos " << globalPos << " has " << (intervalInfo.lastInterval-intervalInfo.firstInterval+1) << " intervals" << endl;
        }
        globalPos += l;
        contigIntervalInfo_.push_back( intervalInfo );
    }
    return intervalCount;
}


// Class RandomIntervalGeneratorUsingIntervalLengthDistribution

pair< unsigned long, unsigned int > RandomIntervalGeneratorUsingIntervalLengthDistribution::getNext( const signed long testValue )
{
    assert( testValue < 0 ); // positive test case not implemented

    if (randomSeq_.hasFinished())
    {
        return make_pair(0,0);
    }
    double randomPos = randomSeq_.getNextAsDouble();

    // randomPos => currentPos
    double probaOfCurrentPos = getIntervalsProbaAtPos( currentPos_ );
    while (randomPos >= accumulatedProbasUntilCurrentPos_ + probaOfCurrentPos)
    {
        accumulatedProbasUntilCurrentPos_ += probaOfCurrentPos;
        ++currentPos_;
        probaOfCurrentPos = getIntervalsProbaAtPos( currentPos_ );
    }

    // randomPos(-currentPos) => fragmentLength
    unsigned int fragmentLength = fragmentLengthDist_.min();
    double probaDiff = randomPos - accumulatedProbasUntilCurrentPos_;
    while (probaDiff > fragmentLengthDist_[fragmentLength])
    {
        probaDiff -= fragmentLengthDist_[fragmentLength++];
    }

    return make_pair( currentPos_, fragmentLength );
}

double RandomIntervalGeneratorUsingIntervalLengthDistribution::getIntervalsProbaAtPos( unsigned long globalPos )
{
    assert( contigLengths_.size() == 1 ); // We restrict the current implementation to one contig at a time
    const unsigned long contigLength = contigLengths_[0];

    if (globalPos + fragmentLengthDist_.max() < contigLength)
    {
        return 1.0;
    }
    else
    {
        unsigned int distanceToContigEnd = contigLength - globalPos;
        unsigned int index = fragmentLengthDist_.max() - distanceToContigEnd;
        if (index >= probas_.size())
        {
            assert( index == probas_.size() ); // We expect to build the index 1 entry at a time
            double p = 0;
            for (unsigned int i=fragmentLengthDist_.min(); i<=distanceToContigEnd; ++i)
            {
                p += fragmentLengthDist_[i];
            }
            probas_.push_back( p );
//            clog << (boost::format("Adding proba %f at index %d (for globalPos %d)") % p % index % globalPos).str() << endl;
        }
        return probas_[index];
    }
}

double RandomIntervalGeneratorUsingIntervalLengthDistribution::getTotalIntervalsProba()
{
    assert( contigLengths_.size() == 1 ); // We restrict the current implementation to one contig at a time
    const unsigned long contigLength = contigLengths_[0];

    double intervalCount = 0;
    for (unsigned long globalPos=0; globalPos<contigLength; ++globalPos)
    {
        intervalCount += getIntervalsProbaAtPos( globalPos );
    }
    if (verbose_)
    {
        EAGLE_DEBUG( 0, "Contig's total intervals probability: " << intervalCount );
    }
    return intervalCount;
}


// Class UniformIntervalGenerator

UniformIntervalGenerator::UniformIntervalGenerator( const vector<unsigned long>& contigLengths, const unsigned int medianFragmentLength, const double step, unsigned long& readCount, bool verbose)
    : IntervalGenerator( contigLengths, readCount, medianFragmentLength, medianFragmentLength, medianFragmentLength, verbose)
    , step_( step )
{
    unsigned long intervalCount = 0;
    unsigned long globalPos = 0;
    BOOST_FOREACH( const unsigned long l, contigLengths_ )
    {
        struct ContigIntervalInfo intervalInfo;
        intervalInfo.firstGlobalPos = globalPos;

        unsigned long contigValidPosCount = l - medianFragmentLength + 1;
        intervalCount += contigValidPosCount;
        intervalInfo.lastGlobalPos = globalPos + contigValidPosCount - 1;

        globalPos += l;
        contigIntervalInfo_.push_back( intervalInfo );
    }

    currentContigIntervalInfo_ = contigIntervalInfo_.begin();
    currentGlobalPos_ = 0 - step_;

    unsigned long newReadCount = static_cast<unsigned long>( (double)intervalCount / step_ );
    if (verbose_)
    {
        clog << (boost::format("Uniform coverage attempts to achieve the specified coverage depth for all the chromosome positions, which is impossible to achieve at the chromosome extremities => the average coverage depth will be lower than specified. Changing read count from %d to %d, for a step of %d") % readCount % newReadCount % step_).str() << endl;
    }
    readCount = newReadCount;
}

pair< unsigned long, unsigned int > UniformIntervalGenerator::getNext( const signed long testValue )
{
    currentGlobalPos_ += step_;
    // Convert intervalNum to interval
    while ( floor(currentGlobalPos_) > currentContigIntervalInfo_->lastGlobalPos )
    {
        if ( currentContigIntervalInfo_+1 != contigIntervalInfo_.end() )
        {
            ++currentContigIntervalInfo_;
            currentGlobalPos_ = currentContigIntervalInfo_->firstGlobalPos;
            if (verbose_)
            {
                clog << "Generating fragments for chromosome starting at global pos " << currentContigIntervalInfo_->firstGlobalPos << endl;
            }
        }
        else
        {
            return make_pair( 0, 0 );
            //EAGLE_WARNING( "Too many reads generated at the end!" );
        }
    }
    assert( floor(currentGlobalPos_) >= currentContigIntervalInfo_->firstGlobalPos );

    return std::make_pair( (unsigned long)floor(currentGlobalPos_), medianFragmentLength_ );
}


pair< unsigned long, unsigned int > RandomIntervalGeneratorFromProbabilityMatrix::getNext( const signed long testValue )
{
    if (remainingSampleCount_ > 0)
    {
        double rand0to1 = (1.0+rand())/(1.0+RAND_MAX);
        curMax_ *= exp(log(rand0to1)/remainingSampleCount_);
        remainingSampleCount_--;
        double choice = (1-curMax_) * intervalSize_;

        lastChoiceRemainder += choice - lastChoice;
        while (lastChoiceRemainder >= fragmentLengthProbabilityMatrix_.getProbabilities()[ lastChoiceIndex ])
        {
            lastChoiceRemainder -= fragmentLengthProbabilityMatrix_.getProbabilities()[ lastChoiceIndex ];
            if (lastChoiceIndex < fragmentLengthProbabilityMatrix_.getProbabilities().size())
            {
                ++lastChoiceIndex;
            }
            else
            {
                // Avoid going past the end of the array because of rounding errors
                break;
            }
        }
        lastChoice = choice;
        unsigned long pickedGlobalPos = lastChoiceIndex / fragmentLengthDist_size;
        unsigned long pickedFragmentLength = lastChoiceIndex % fragmentLengthDist_size + fragmentLengthDist_min;
        cout << (boost::format("remainingSampleCount_=%d, choice=%f, choiceIndex=%d, pickedGlobalPos=%d, pickedFragmentLength=%d") % remainingSampleCount_ % choice % lastChoiceIndex % pickedGlobalPos % pickedFragmentLength).str() << endl;
        return make_pair( pickedGlobalPos, pickedFragmentLength );
    }
    else
    {
        EAGLE_WARNING( "Too many reads generated at the end!" );
        return make_pair( 0, 0 );
    }
}


} // namespace model
} // namespace eagle
