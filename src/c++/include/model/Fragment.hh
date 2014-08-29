/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Set of classes to represent DNA fragments
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MODEL_FRAGMENT_HH
#define EAGLE_MODEL_FRAGMENT_HH

#include <utility>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>


using namespace std;


namespace eagle
{
namespace model
{


class Fragment
{
public:
    Fragment(unsigned long startPos, unsigned long fragmentLength, unsigned long fragmentNum) : startPos_(startPos), fragmentLength_(fragmentLength), fragmentNum_(fragmentNum), multiplexedDatasetId(0) {}

    // Useful for vector<Fragment>
    Fragment() : startPos_(0), fragmentLength_(0), fragmentNum_(0), multiplexedDatasetId(0) {}
    bool isValid() { return (fragmentLength_ > 0); }
    void invalidate() { fragmentLength_ = 0; }

    friend std::ostream& operator<<( std::ostream& os, const Fragment& f )
    {
        return os << (boost::format("{%d,%d}") % f.startPos_ % f.fragmentLength_).str();
    }

    double getGcContent();

    unsigned long startPos_;
    unsigned long fragmentLength_;
    unsigned long fragmentNum_;
    unsigned int multiplexedDatasetId;    //    MultiplexedDatasetInfo &provenance; // TODO: refine this
};


class FragmentWithAllocationMetadata : public Fragment
{
public:
    FragmentWithAllocationMetadata();
    FragmentWithAllocationMetadata( const std::pair<unsigned long,unsigned int>& p );
    FragmentWithAllocationMetadata(unsigned long startPos, unsigned long fragmentLength, unsigned long fragmentNum, unsigned int tile) : Fragment( startPos, fragmentLength, fragmentNum ), allocatedTile_( tile ) {}
    void allocateRandomTile(unsigned long tileCount);
    void allocateTileInSequence(unsigned long tileCount, unsigned long readNum, unsigned long readCount);
    void allocateInterleavedTile(unsigned long tileCount);
    bool operator<( const FragmentWithAllocationMetadata& rhs ) const;

    friend std::ostream& operator<<( std::ostream& os, const FragmentWithAllocationMetadata& f )
    {
        return os << (boost::format("{%d,%d,%d}") % f.startPos_ % f.fragmentLength_ % f.allocatedTile_).str();
    }

    unsigned int allocatedTile_;
};


class FragmentList
{
public:
    FragmentList( const boost::filesystem::path& dir = "", const unsigned long firstRequestedPos = 0, const unsigned long lastRequestedPos = 0xFFFFFFFFFFFFFFFFul, const unsigned long fetchBefore = 0 );
    Fragment getNext( unsigned int tilePattern, unsigned int mask, unsigned int *tilePtr=0 ); // from file
    inline Fragment getNext( unsigned int desiredTile ) { return getNext( desiredTile, -1 ); } // Get fragments from specific tile
    inline Fragment getNext() { return getNext( 0, 0 ); } // Get fragments from any tile
    inline bool getNext( Fragment& resultFragment ) { resultFragment = getNext(); return resultFragment.isValid(); }
    FragmentWithAllocationMetadata getNextWithTile( unsigned int tilePattern = 0, unsigned int mask = -1, unsigned int *tilePtr=0 );
    unsigned long long getTileSize( unsigned int tileNum );
    unsigned long long size();
    
private:
    ifstream in1, in2, in3, in4;
    unsigned int fragmentNum_;


    unsigned long insertSize_;

    unsigned int currentContigNum_;
    unsigned long currentContigLength_;
    unsigned long currentPos_;

    unsigned long firstRequestedPos_;
    unsigned long lastRequestedPos_;
};


class MultiFragmentFilesReader
{
public:
    MultiFragmentFilesReader( const std::vector<unsigned long>& contigLengths, const std::vector<string>& contigNames, const boost::filesystem::path& dir, unsigned long& readCount, const bool verbose=true );
    virtual ~MultiFragmentFilesReader() {}
    FragmentWithAllocationMetadata getNext( const signed long testValue = -1 );

private:
    bool openNextFragmentList( bool advance );
/*
    double step_;
    struct ContigIntervalInfo
    {
        unsigned long firstGlobalPos;
        unsigned long lastGlobalPos;
    };
    std::vector<ContigIntervalInfo> contigIntervalInfo_;
    std::vector<ContigIntervalInfo>::iterator currentContigIntervalInfo_;
*/
    const std::vector<unsigned long> contigLengths_;
    const std::vector<string> contigNames_;
    const boost::filesystem::path dir_;
    unsigned int currentContigNum_;
    unsigned long currentContigStartPos_;
    FragmentList *fragmentList_;
};



} // namespace model
} // namespace eagle

#endif // EAGLE_MODEL_FRAGMENT_HH
