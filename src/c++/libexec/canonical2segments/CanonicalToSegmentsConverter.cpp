/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include <boost/ptr_container/ptr_deque.hpp>
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "genome/VariantList.hh"
#include "genome/ReferenceToSample.hh"
#include "CanonicalToSegmentsConverter.hh"


using namespace std;


namespace eagle
{
namespace main
{

struct AlleleEvent
{
    AlleleEvent(
        string& chrAllele, 
        unsigned long pos,
        genome::Event* event
        )
        : chrAllele_( chrAllele )
        , pos_( pos )
        , event_( event )
        {
        }

    bool operator<( const AlleleEvent& rhs ) const
    {
        int comp = chrAllele_.compare( rhs.chrAllele_ );
        if (comp == 0)
        {
            return (pos_ < rhs.pos_);
        }
        else
        {
            return comp < 0;
        }
    }

    friend std::ostream& operator<<( std::ostream& os, const AlleleEvent& obj )
    {
        return os << "{ " << obj.chrAllele_ << ", " << obj.pos_ << ", " << *(obj.event_) << " } ";
    }

    string chrAllele_;
    unsigned long pos_;
    genome::Event* event_;
};




CanonicalToSegmentsConverter::CanonicalToSegmentsConverter (const CanonicalToSegmentsConverterOptions &options )
    : options_( options )
{
}

/*
class LtCompareByDestInfo
{
public:
    LtCompareByDestInfo( const string& chrAllele )
    {
    }

    bool operator()( const eagle::genome::Event& e1, const eagle::genome::Event &e2 )
    {
        if( e1.metadata_.hasInfo() && event->metadata_.getInfo("DEST").size() > 0)
            return true;
    }
};
*/

void CanonicalToSegmentsConverter::run()
{
    const vector< boost::filesystem::path > inputs( 1, options_.input );
    const boost::filesystem::path output;
    const model::Ploidy ploidy(1);
    genome::VariantList variantList( inputs, output, ploidy);
//    unsigned long timeProcessing = 0;
    unsigned long timeIO = 0;

    clock_t start = clock();

    clog << "Loading " << variantList.fileCount() << " variant list" << (variantList.fileCount() - 1 ? "s " : " ") << "..." << endl;
    const bool filterSnpOut = true;
    const bool filterBeginEndMarkersOut = false;
    variantList.load( filterSnpOut, filterBeginEndMarkersOut );
    clog << "Loaded " << variantList.size() << " event" << (variantList.size() - 1 ? "s " : " ")
              << "in " << eagle::common::displayTime(clock() - start, timeIO) << endl;

    // Go through the list of variants, and set/mark each variant that has a DEST field as a segment start.
    // Delete the variants that don't have a DEST field
    // Duplicate the variants that have multiple DEST values
    // Output this list of segment start
    vector< unsigned long > breakends;
    vector< AlleleEvent > alleleEvents;
    for (eagle::genome::EventIterator event = variantList.begin(); event != variantList.end(); ++event)
    {
        if( event->metadata_.hasInfo() && event->metadata_.getInfo("DEST").size() > 0)
        {
//            cout << "New segment start: " << *event << endl;
            breakends.push_back( event->from() );
            breakends.push_back( event->to() );
            BOOST_FOREACH( const string& destField, event->metadata_.getInfo("DEST") )
            {
//                cout << destField << endl;
                vector<string> items;
                boost::split( items, destField, boost::is_any_of(":") );
                assert( items.size() == 2 );

                struct AlleleEvent alleleEvent( 
                    items[0]
                    , boost::lexical_cast<unsigned long, string>( items[1] )
                    , &*event
                    );
                alleleEvents.push_back( alleleEvent );
            }
        }
        else
        {
//            cout << "Skipped: " << *event << endl;
/*
            --it;
            variantList.remove( it+1 );
*/
        }
    }
    sort( breakends.begin(), breakends.end() );

    // Sort the variant list by DEST=chromosome_allele:#
    // Go through the sorted list, and for each chromosome allele, extract the list of segments it went through:
    //   For each variant, the segment number is obtained from the variant position, and the list of consecutive segments is calculated by using the distant to the next variant.
    sort( alleleEvents.begin(), alleleEvents.end() );
    vector< genome::RefToSampleSegment > refToSampleSegments;
    string lastSampleChrAllele = "";

    string lastRefChr = "";
    unsigned long lastRefPos = 0;
    unsigned long lastSamplePos = 0;
    BOOST_FOREACH( const AlleleEvent& alleleEvent, alleleEvents )
    {
        if (alleleEvent.chrAllele_ != lastSampleChrAllele)
        {
            cout << "New sample chr allele: " << alleleEvent.chrAllele_  << endl;
            lastSampleChrAllele = alleleEvent.chrAllele_;
            lastSamplePos = 0;
            lastRefChr = alleleEvent.event_->src();
            // TODO: the following line is wrong if the Genome Mutator's starting point was at the end of a ref chromosome
            lastRefPos = 0;
        }
        cout << alleleEvent << endl;

        genome::RefToSampleSegment refToSampleSegment;

        unsigned long endSamplePos = alleleEvent.pos_;
        unsigned long newSamplePos = alleleEvent.pos_ +  alleleEvent.event_->getVariant().sequence.size();
        signed long sampleFragmentLength = endSamplePos - lastSamplePos;
        refToSampleSegment.sampleChrAllele_ = lastSampleChrAllele;
        refToSampleSegment.samplePos_ = lastSamplePos + 1;
        cout << (boost::format("  Sample: %s:%d->%d (%d)") % lastSampleChrAllele % lastSamplePos % endSamplePos % sampleFragmentLength).str() << endl;
        lastSamplePos = newSamplePos;

        unsigned long endRefPos = alleleEvent.event_->to();
        string newRefChr = alleleEvent.event_->dest();
        unsigned long newRefPos = alleleEvent.event_->from();
        signed long refFragmentLength = endRefPos - lastRefPos;
        refToSampleSegment.refChr_ = lastRefChr;
        refToSampleSegment.refPos_ = min( endRefPos, lastRefPos ) + 1;
        refToSampleSegment.segmentLengthWithRefDirection_ = refFragmentLength;
        cout << (boost::format("  Ref   : %s:%d->%d (%d) - new=%s:%d") % lastRefChr % lastRefPos % endRefPos % refFragmentLength % newRefChr % newRefPos).str() << endl;
        lastRefChr = newRefChr;
        lastRefPos = newRefPos;

        assert( sampleFragmentLength >= 0 );
        assert( sampleFragmentLength == labs(refFragmentLength) || sampleFragmentLength == labs(refFragmentLength)-1 );

        if( sampleFragmentLength > 0 )
        {
            if (alleleEvent.event_->isBeginEndMarker())
            {
                // The end marker of a chromosome is 1 base after the end of the chromosome => the length of the segment needs to be adjusted
                --(refToSampleSegment.segmentLengthWithRefDirection_);
            }
            refToSampleSegments.push_back( refToSampleSegment );
        }
    }


    // Sort the list of RefToSampleSegment segments by ref position, and output them
    stable_sort( refToSampleSegments.begin(), refToSampleSegments.end() ); // Using stable_sort to keep allele1 before allele2 when they happen at the same ref position
    {
        boost::filesystem::path outFilename;
        if (options_.outputDir != "")
        {
            outFilename = options_.outputDir / "segmentsFromRef.tsv";
        }
        else
        {
            outFilename = "segmentsFromRef.tsv";
        }
        clog << "Output: " << outFilename.string() << endl;
        ofstream refToSampleSegmentOutFile( outFilename.string().c_str() );
        assert( refToSampleSegmentOutFile.good() && "Error opening the output file in writing" );

        refToSampleSegmentOutFile << "refChr\tleftMostRefPos\tsampleChrAllele\tsamplePos\tlengthWithRefDirection" << endl;
        BOOST_FOREACH( const genome::RefToSampleSegment& refToSampleSegment, refToSampleSegments )
        {
            refToSampleSegmentOutFile << refToSampleSegment << endl;
        }
    }
}


} // namespace main
} // namespace eagle
