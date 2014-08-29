/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Set of classes to manipulate conversions between reference and sample genome coordinates 
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_GENOME_REFERENCE_TO_SAMPLE_HH
#define EAGLE_GENOME_REFERENCE_TO_SAMPLE_HH

#include <utility>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "genome/SharedFastaReference.hh"

using namespace std;


namespace eagle
{
namespace genome
{


struct RefToSampleSegment
{
    string refChr_;
    unsigned long refPos_;

    string sampleChrAllele_;
    unsigned long samplePos_;

    long segmentLengthWithRefDirection_;

    RefToSampleSegment* previousSegmentInGroup_;
    RefToSampleSegment* nextSegmentInGroup_;

    RefToSampleSegment()
        : refPos_(0)
        , samplePos_(0)
        , segmentLengthWithRefDirection_(0)
        , previousSegmentInGroup_(0)
        , nextSegmentInGroup_(0)
        {}

    bool operator<( const RefToSampleSegment& rhs ) const
    {
        int comp = refChr_.compare( rhs.refChr_ );
        if (comp == 0)
        {
            return (refPos_ < rhs.refPos_);
        }
        else
        {
            return comp < 0;
        }
    }

    unsigned long getRightMostRefPos()
    {
        return refPos_ + labs(segmentLengthWithRefDirection_) - 1;
    }

    unsigned long getSampleGlobalStartPos()
    {
        eagle::model::Locus localPos( sampleChrAllele_, samplePos_ );
        unsigned long globalPos = genome::SharedFastaReference::get()->local2global( localPos );
        return globalPos;
    }

    unsigned long getSampleGlobalEndPos()
    {
        eagle::model::Locus localPos( sampleChrAllele_, samplePos_ + labs(segmentLengthWithRefDirection_) - 1 );
        unsigned long globalPos = genome::SharedFastaReference::get()->local2global( localPos );
        return globalPos;
    }

    friend std::ostream& operator<<( std::ostream& os, const RefToSampleSegment& obj )
    {
        return os << (boost::format("%s\t%d\t%s\t%d\t%d")
                      % obj.refChr_
                      % obj.refPos_
                      % obj.sampleChrAllele_
                      % obj.samplePos_
                      % obj.segmentLengthWithRefDirection_
            ).str();
    }
};

/*
class RefToSampleSegmentGroup : public RefToSampleSegment
{
public:
    RefToSampleSegmentGroup( ) {}

    RefToSampleSegmentGroup( const RefToSampleSegment& firstSegment ) : RefToSampleSegment( firstSegment )
    {
        lastAddedSegment = firstSegment;
    }

    void tryToAdd( const RefToSampleSegment& segment , bool& needsDifferentGroup, bool& needsNewGroup )
    {
        needsDifferentGroup = needsNewGroup = false;
        if (segment.sampleChrAllele_ != sampleChrAllele_)
        {
            needsDifferentGroup = true;
        }
        else
        {
            bool isIns = (lastAddedSegment.samplePos_ + lastAddedSegment.segmentLengthWithRefDirection_ != segment.samplePos_);
            bool isDel = (lastAddedSegment.refPos_ + lastAddedSegment.segmentLengthWithRefDirection_ != segment.refPos_);
            if (isIns && isDel)
            {
                needsNewGroup = true;
            }
            else
            {
                segmentLengthWithRefDirection_ += addEntryToCigarModifier( lastAddedSegment, segment );
                segmentLengthWithRefDirection_ += segment.segmentLengthWithRefDirection_;
                lastAddedSegment = segment;
            }
        }
    }

    std::string getCigarModifier()
    {
        return cigarModifier_;
    }

private:
    int addEntryToCigarModifier( const RefToSampleSegment& s1, const RefToSampleSegment& s2 )
    {
        int insLength = (s1.samplePos_ + s1.segmentLengthWithRefDirection_ - s2.samplePos_);
        int delLength = (s1.refPos_ + s1.segmentLengthWithRefDirection_ - s2.refPos_);
        if (insLength && delLength)
        {
            assert( false && "should never reach here" );
        }
        else
        {
            if (insLength)
            {
                cigarModifier_ += (boost::format("%dM%dI") % s1.segmentLengthWithRefDirection_ % insLength).str();
            }
            else if (delLength)
            {
                cigarModifier_ += (boost::format("%dM%dD") % s1.segmentLengthWithRefDirection_ % delLength).str();
            }
            else
            {
                // end of chromosome?
                assert( false && "should never reach here - strange segment transition without indel" );
            }
        }
        return 0;
    }

    RefToSampleSegment lastAddedSegment;
    string cigarModifier_;
};
*/

class RefToSampleSegmentReader
{
public:
/*
    RefToSampleSegmentReader( const boost::filesystem::path& filename ) : ifs_(filename.string().c_str()), isCacheUsed_(false)
    {
        clog << "+ Input segment map: \"" << filename << "\"" << endl;
        if (!ifs_.good())
        {
            EAGLE_ERROR( (boost::format("Error opening file %s") % filename.string()).str() );
        }
    }
*/
    RefToSampleSegmentReader( const boost::filesystem::path& filename, const string& requestedChr ) : ifs_(filename.string().c_str()), isCacheUsed_(false), index_(0)
    {
        clog << "+ Input segment map: \"" << filename << "\"" << endl;
        if (!ifs_.good())
        {
            EAGLE_ERROR( (boost::format("Error opening file %s") % filename.string()).str() );
        }

        // Read all the segments
        RefToSampleSegment refToSampleSegment;
        while (getNextSegmentForRefChr( requestedChr, refToSampleSegment ))
        {
            segments_.push_back( refToSampleSegment );
        }

        // Identify segment groups
        std::vector<bool> segmentConsumed( segments_.size(), false );
        for (unsigned int i = 0; i<segments_.size(); ++i)
        {
            if (!segmentConsumed[i])
            {
                // create a new group starting from this segment
//                RefToSampleSegmentG_____ newGroup( segments_[i] );
                int previousSegmentIndex = i;
                segmentConsumed[i] = true;
                unsigned int savedIndex = i;
                for (++i; i<segments_.size(); ++i)
                {
                    if (!segmentConsumed[i])
                    {
                        bool needsDifferentGroup, needsNewGroup;
                        checkIfSameGroup( segments_[previousSegmentIndex], segments_[i], needsDifferentGroup, needsNewGroup );
                        if (needsNewGroup)
                        {
                            break; // the new group will be created at the next iteration, which we try to start just after the current group begins
                        }
                        else if (needsDifferentGroup)
                        {
                            // do nothing, the appropriate group will be created in another iteration
                        }
                        else
                        {
                            segments_[previousSegmentIndex].nextSegmentInGroup_ = &segments_[i];
                            segments_[i].previousSegmentInGroup_ = &segments_[previousSegmentIndex];
                            segmentConsumed[i] = true;
                            previousSegmentIndex = i;
                        }
                    }
                }
                // flush this group out
//                segmentGroups_.push_back( newGroup );
                i = savedIndex;
            }
        }
    }

    void checkIfSameGroup( const RefToSampleSegment& s1, const RefToSampleSegment& s2, bool& needsDifferentGroup, bool& needsNewGroup )
    {
        needsDifferentGroup = needsNewGroup = false;
        if (s1.sampleChrAllele_ != s2.sampleChrAllele_)
        {
            needsDifferentGroup = true;
        }
        else
        {
            bool isIns = (s1.samplePos_ + s1.segmentLengthWithRefDirection_ != s2.samplePos_);
            bool isDel = (s1.refPos_ + s1.segmentLengthWithRefDirection_ != s2.refPos_);
            if (isIns && isDel)
            {
//                needsNewGroup = true;
            }
        }
    }

/*
    bool getNextSegmentGroup( RefToSampleSegmentGroup& refToSampleSegment )
    {
        if (groupIndex_ < segmentGroups_.size())
        {
             refToSampleSegment = segmentGroups_[groupIndex_++];
             return true;
        }
        return false;
    }
*/
/*
    bool lookAheadNextSegment( RefToSampleSegment& refToSampleSegment )
    {
        
    }
*/
    bool getNextSegment( RefToSampleSegment& refToSampleSegment )
    {
        if (index_ < segments_.size())
        {
             refToSampleSegment = segments_[index_++];
             return true;
        }
        return false;
    }
    void goBack( const unsigned int distance )
    {
        assert( index_ >= distance );
        index_ -= distance;
    }

    bool getNextSegmentForRefChr( const string& requestedChr, RefToSampleSegment& refToSampleSegment )
    {
        if (isCacheUsed_)
        {
            refToSampleSegment = cachedSegment_;
            isCacheUsed_ = false;
            return true;
        }
        else
        {
            string line;
            while (true)
            {
                std::getline( ifs_, line );
                if (ifs_.good())
                {
                    vector<string> items;
                    boost::split( items, line, boost::is_any_of("\t") );
                    assert( items.size() == 5 );
                    cachedSegment_.refChr_ = items[0];
                    if (cachedSegment_.refChr_ != requestedChr)
                    {
                        continue;
                    }
                    cachedSegment_.refPos_ = boost::lexical_cast<unsigned long, string>( items[1] );
                    cachedSegment_.sampleChrAllele_ = items[2];
                    cachedSegment_.samplePos_ = boost::lexical_cast<unsigned long, string>( items[3] );
                    cachedSegment_.segmentLengthWithRefDirection_ = boost::lexical_cast<signed long, string>( items[4] );
                    refToSampleSegment = cachedSegment_;
                    // Note: isCachedUsed==false at this stage. We do not set it until the user calls the goBack() method
                    return true;
                }
                else
                {
                    return false;
                }
            }
        }
    }
/*
    void goBack( const unsigned int distance )
    {
        assert( distance == 1 );
        assert( !isCacheUsed_ );
        isCacheUsed_ = true;
    }
*/

private:
    ifstream ifs_;
    RefToSampleSegment cachedSegment_;
    bool isCacheUsed_;

    std::vector<RefToSampleSegment> segments_;
    unsigned int index_;
//    std::vector<RefToSampleSegmentGroup> segmentGroups_;
//    unsigned int groupIndex_;
};


} // namespace genome
} // namespace eagle

#endif // EAGLE_GENOME_REFERENCE_TO_SAMPLE_HH
