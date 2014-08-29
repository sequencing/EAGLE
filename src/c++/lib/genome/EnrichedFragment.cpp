/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "genome/SharedFastaReference.hh"
#include "genome/Reference.hh"
#include "genome/EnrichedFragment.hh"

using namespace std;


namespace eagle
{
namespace genome
{


FragmentComponentHardcoded::FragmentComponentHardcoded( string bases)
    : FragmentComponent()
    , bases_( bases )
{
    length_ = bases_.size();
}

char FragmentComponentHardcoded::getBase( const unsigned int posInRead, const eagle::model::Fragment &fragment, const bool isForward ) const
{
    //        const unsigned long fragmentStartPos = fragment.startPos_;
    char base;
    assert( posInRead < bases_.size() );
    base = bases_[posInRead];
    
    if (!isForward)
    {
        // reverse complement
        base = ~base;
        //bclBase = (bclBase & ~3) | (~bclBase & 3); 
    }
    
    return base;
}


FragmentComponentRealDna::FragmentComponentRealDna()
    : FragmentComponent()
{
    length_ = 0;
}

char FragmentComponentRealDna::getBase( const unsigned int posInRead, const eagle::model::Fragment &fragment, const bool isForward ) const
{
    const unsigned long fragmentStartPos = fragment.startPos_;
    const unsigned long fragmentLength = fragment.fragmentLength_;
    char base;

    bool overlapContigBoundary = false;
    if (isForward)
    {
        base = eagle::genome::SharedFastaReference::get()->get( fragmentStartPos, posInRead, overlapContigBoundary );
    }
    else
    {
        base = eagle::genome::SharedFastaReference::get()->get( fragmentStartPos, fragmentLength - static_cast<unsigned long>(posInRead) - 1, overlapContigBoundary );
        // reverse complement
        base = ~base;
    }
    if (overlapContigBoundary)
    {
        EAGLE_ERROR( (boost::format( "overlapContigBoundary == true : fragmentStartPos=%d, fragmentLength=%d, posInRead=%d, isForward=%d (TODO) This assert usually reveals a bug, but is sometimes due to the not-yet-implemented features to continue reads on adapters when the DNA ends") % fragmentStartPos % fragmentLength % posInRead % isForward).str() );
    }

    return base;
}


char FragmentStructure::getBase( const unsigned int readNum, const unsigned int posInRead, const eagle::model::Fragment &fragment ) const
{
    assert( readNum < reads_.size() );
    unsigned int compNum = reads_[readNum].first;
    bool isForward = reads_[readNum].second;
    char base = components_[compNum]->getBase( posInRead, fragment, isForward );
    return base;
}

unsigned int FragmentStructure::getReadLength( const unsigned int readNum, const eagle::model::Fragment &fragment ) const
{
    assert( readNum < reads_.size() );
    unsigned int compNum = reads_[readNum].first;
    unsigned int result = components_[compNum]->length_;
    if (result == 0)
    {
        result = fragment.fragmentLength_;
    }
    return result;
}


EnrichedFragment::EnrichedFragment( const eagle::model::Fragment &fragment, const vector<FragmentStructure> &multiplexedFragmentStructures, const unsigned int dnaFragmentDirection )
    : fragment_( fragment )
    , dnaFragmentDirection_( dnaFragmentDirection )
    , structure_( multiplexedFragmentStructures[fragment.multiplexedDatasetId*2 + dnaFragmentDirection_] )
{
}

char EnrichedFragment::getBase( unsigned int read, unsigned int posInRead ) const
{
    char base = structure_.getBase( read, posInRead, fragment_ );
    return base;
}

unsigned int EnrichedFragment::getReadCount() const
{
    return structure_.reads_.size();
}

unsigned int EnrichedFragment::getReadLength( unsigned int r ) const
{
    return structure_.getReadLength( r, fragment_ ); 
}


} // namespace genome
} // namespace eagle
