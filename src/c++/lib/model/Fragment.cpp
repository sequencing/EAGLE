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
#include "model/Fragment.hh"


using namespace std;


namespace eagle
{
namespace model
{

/*
double Fragment::getGcContent()
{
    EAGLE_WARNING( "TODO: Fragment::getGcContent()" );
    return 0.5;
}
*/

FragmentWithAllocationMetadata::FragmentWithAllocationMetadata()
    : Fragment(0,0,0)
    , allocatedTile_(0)
{
}

FragmentWithAllocationMetadata::FragmentWithAllocationMetadata( const pair<unsigned long,unsigned int>& p )
    : Fragment( p.first, p.second, 0 )
    , allocatedTile_( 0 )
{
}

void FragmentWithAllocationMetadata::allocateRandomTile(unsigned long tileCount)
{
    allocatedTile_ = rand() % tileCount;
}

void FragmentWithAllocationMetadata::allocateTileInSequence(unsigned long tileCount, unsigned long readNum, unsigned long readCount)
{
    allocatedTile_ = (tileCount * readNum) / readCount;
}

void FragmentWithAllocationMetadata::allocateInterleavedTile(unsigned long tileCount)
{
    static unsigned long nextTile = 0;
    allocatedTile_ = nextTile;
    nextTile = (nextTile+1) % tileCount;
}

bool FragmentWithAllocationMetadata::operator<( const FragmentWithAllocationMetadata& rhs ) const
{
    if (startPos_ < rhs.startPos_) return true;
    if (startPos_ > rhs.startPos_) return false;
    if (fragmentLength_ < rhs.fragmentLength_) return true;
    return false;
}


FragmentList::FragmentList( const boost::filesystem::path& dir, const unsigned long firstRequestedPos, const unsigned long lastRequestedPos, const unsigned long fetchBefore )
    : in1( (dir/"fragments.pos"   ).string().c_str(), ios::binary )
    , in2( (dir/"fragments.length").string().c_str(), ios::binary )
    , in3( (dir/"fragments.tile"  ).string().c_str(), ios::binary )
    , in4( (dir/"fragments.stats" ).string().c_str(), ios::binary )
    , fragmentNum_( 0 )
    , currentPos_( 0 )
    , firstRequestedPos_( firstRequestedPos )
    , lastRequestedPos_( lastRequestedPos )
{
    unsigned long startPos = (firstRequestedPos > fetchBefore) ? (firstRequestedPos - fetchBefore) : 0;

    // If a non-zero position is requested, use index file to jump as close as possible
    if (startPos > 0)
    {
        // Reading fragments.pos.index v1
        ifstream indexFile( (dir/"fragments.pos.index").string().c_str() );
        unsigned long version;
        unsigned long indexInterval;
        indexFile.read( (char*)&version, sizeof(unsigned long));
        assert( version == 1 );
        indexFile.read( (char*)&indexInterval, sizeof(unsigned long));
        unsigned long pos = 0;
        unsigned long previousPos = 0;
        while (indexFile.good() && pos < startPos)
        {
            previousPos = pos;
            indexFile.read( (char*)&pos, sizeof(unsigned long));
        }
        if (!indexFile.good())
        {
            // if we reached the end, move properly to the end, so that tellg returns the file size rather than -1
            indexFile.clear();
            indexFile.seekg( 0, ios_base::end );
        }
        unsigned long posInIndexFile = indexFile.tellg();
        unsigned long indexEntryNum = posInIndexFile / sizeof(unsigned long) - 3;
        fragmentNum_ = indexEntryNum * indexInterval;

        ifstream shiftFile( (dir/"fragments.pos.shift").string().c_str() );
        unsigned int shift;
        shiftFile.seekg( indexEntryNum * sizeof(unsigned int) );
        shiftFile.read( (char*)&shift, sizeof(unsigned int));

        unsigned long posInPosFile =  (fragmentNum_ + shift) * 2;
        unsigned long posInOtherFiles =  fragmentNum_ * 2;
        in1.seekg( posInPosFile );
        in2.seekg( posInOtherFiles );
        in3.seekg( posInOtherFiles );
        currentPos_ = previousPos;
    }
}

Fragment FragmentList::getNext( unsigned int desiredTile, unsigned int mask, unsigned int *tilePtr )
{
    unsigned long length=0;
    unsigned int tile=0;
    do {
        /* Format: binary 6 bytes per fragment */
        unsigned long long posDiff = 0;
        in1.read( (char*)&posDiff, 2);
        in2.read( (char*)&length, 2);
        in3.read( (char*)&tile, 2);
        if (!in1.good())
        {
            return Fragment(); // returns a fragment that has .isValid()==false
        }
        assert( in2.good() );
        assert( in3.good() );
        if (posDiff == 65535)
        {
            unsigned long posDiffByte2 = 0;
            unsigned long posDiffByte1 = 0;
            unsigned long posDiffByte0 = 0;
            in1.read( (char*)&posDiffByte2, 2);
            in1.read( (char*)&posDiffByte1, 2);
            in1.read( (char*)&posDiffByte0, 2);
            posDiff = (posDiffByte2 << 32) | (posDiffByte1 << 16) | posDiffByte0;
        }
        currentPos_ += posDiff;
        ++fragmentNum_;
    } while ( (tile & mask) != desiredTile || ((currentPos_ < firstRequestedPos_) && (currentPos_+length-1 < firstRequestedPos_)) );

    if (currentPos_ > lastRequestedPos_)
    {
        return Fragment(); // returns a fragment that has .isValid()==false
    }

    if (tilePtr)
    {
        *tilePtr = tile;
    }

    Fragment fragment( currentPos_, length, fragmentNum_-1 );
    return fragment;
}

FragmentWithAllocationMetadata FragmentList::getNextWithTile( unsigned int desiredTile, unsigned int mask, unsigned int *tilePtr )
{
    unsigned long length=0;
    unsigned int tile=0;
    do {
        /* Format: binary 6 bytes per fragment */
        unsigned long long posDiff = 0;
        in1.read( (char*)&posDiff, 2);
        in2.read( (char*)&length, 2);
        in3.read( (char*)&tile, 2);
        if (!in1.good())
        {
            return FragmentWithAllocationMetadata(); // returns a fragment that has .isValid()==false
        }
        assert( in2.good() );
        assert( in3.good() );
        if (posDiff == 65535)
        {
            unsigned long posDiffByte2 = 0;
            unsigned long posDiffByte1 = 0;
            unsigned long posDiffByte0 = 0;
            in1.read( (char*)&posDiffByte2, 2);
            in1.read( (char*)&posDiffByte1, 2);
            in1.read( (char*)&posDiffByte0, 2);
            posDiff = (posDiffByte2 << 32) | (posDiffByte1 << 16) | posDiffByte0;
        }
        currentPos_ += posDiff;
        ++fragmentNum_;
    } while ( (tile & mask) != desiredTile || ((currentPos_ < firstRequestedPos_) && (currentPos_+length-1 < firstRequestedPos_)) );

    if (currentPos_ > lastRequestedPos_)
    {
        return FragmentWithAllocationMetadata(); // returns a fragment that has .isValid()==false
    }

    if (tilePtr)
    {
        *tilePtr = tile;
    }

    FragmentWithAllocationMetadata fragment( currentPos_, length, fragmentNum_-1, tile );
    return fragment;
}

unsigned long long FragmentList::getTileSize( unsigned int tileNum )
{
    unsigned int tileReadCount = 0;

    in4.clear();
    in4.seekg( tileNum * sizeof( unsigned int ) );
    in4.read( (char*)&tileReadCount, sizeof( unsigned int ) );
    return tileReadCount;
}

unsigned long long FragmentList::size()
{
    unsigned long long allTilesReadCount = 0;

    in4.clear();
    in4.seekg( 0 );
    while (in4.good())
    {
        unsigned int tileReadCount = 0;
        in4.read( (char*)&tileReadCount, sizeof( unsigned int ) );
        allTilesReadCount += tileReadCount;
    }
    return allTilesReadCount;
}


//class MultiFragmentFilesReader
MultiFragmentFilesReader::MultiFragmentFilesReader( const std::vector<unsigned long>& contigLengths, const std::vector<string>& contigNames, const boost::filesystem::path& dir, unsigned long& readCount, const bool verbose )
//    : IntervalGenerator( contigLengths, readCount, 0, 0, 0, verbose)
    : contigLengths_( contigLengths )
    , contigNames_( contigNames )
    , dir_( dir )
    , currentContigNum_( 0 )
    , currentContigStartPos_( 0 )
    , fragmentList_( 0 )
{
    openNextFragmentList( false );
}

FragmentWithAllocationMetadata MultiFragmentFilesReader::getNext( const signed long testValue )
{
    while (1)
    {
        FragmentWithAllocationMetadata f = fragmentList_->getNextWithTile( 0, 0 );
        if (f.isValid())
        {
            f.startPos_ += currentContigStartPos_;
            return f;
        }
        else
        {
            if (!openNextFragmentList( true ))
            {
                return FragmentWithAllocationMetadata();
            }
        }
    }
}

bool MultiFragmentFilesReader::openNextFragmentList( bool advance )
{
    if (fragmentList_)
    {
        delete fragmentList_;
        fragmentList_ = 0;
    }
    if (currentContigNum_ >= contigNames_.size())
    {
        return false;
    }
    if (advance)
    {
        currentContigStartPos_ += contigLengths_[ currentContigNum_ ];
        ++currentContigNum_;
        if (currentContigNum_ >= contigNames_.size())
        {
            return false;
        }
    }
    const boost::filesystem::path dir( dir_/(string("fragments_")+contigNames_[ currentContigNum_ ]) );
    if (is_directory(dir))
    {
        clog << "Opening directory " << dir << endl;
        fragmentList_ = new FragmentList( dir );
        return true;
    }
    else
    {
        EAGLE_ERROR( (boost::format("Missing directory %s") % dir.string()).str() );
    }
    return false;
}


} // namespace model
} // namespace eagle
