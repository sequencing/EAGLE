/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Set of classes to read fragments.{pos|tile|...}, cache and access the information efficiently
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MODEL_FRAGMENT_POS_RESOLVER_HH
#define EAGLE_MODEL_FRAGMENT_POS_RESOLVER_HH

#include <utility>
#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include "Fragment.hh"


using namespace std;


namespace eagle
{
namespace model
{

class FragmentPosResolver
{
public:
    FragmentPosResolver() {}
    unsigned long GetSimulatedPosInSampleGenome( const unsigned int lane, const unsigned int tile, const unsigned long cluster )
    {
        clog << (boost::format("FragmentPosResolver::GetSimulatedPosInSampleGenome: lane=%d, tile=%d, cluster=%d") % lane % tile % cluster).str() << endl;        
        unsigned int fullTileNum = (lane-1)*32 +  + (((tile/1000)-1)%10)*16 + (((tile/100)-1)%10)*8 + (tile-1)%10;
        clog << (boost::format("  fullTileNum=%d=0x%x") % fullTileNum % fullTileNum).str() << endl;        

/*
  keep a cache of lane+tile+cluster => global sample pos
    for each non-cached entry: start from closest smaller entry and go through the fragments.tile and fragments.pos files. Update cache

idea 1 (slowest):
 read&accumulate fragment.pos values from the start + count in fragments.tile each time we see our tile number
  until we reach the requested tileFragmentNumber
idea 2 (slow):
 count in fragments.tile each time we see our tile number
  until we reach the requested tileFragmentNumber
 that gives us the globalFragmentNumber
 use fragment.pos.index to find a starting point close to globalFragmentNumber
 finish with read&accumulate fragment.pos from this starting point

idea 3:
 same as idea 2, but avoid counting fragments.tile from the start:
  we need, for each tile, an index "fragmentNumber=>globalFragmentNumber"
   globalFragmentNumber should probably be a subset of those in fragment.pos.index



check v1: (notation: a-%b = a-a%b)
 fragments.tile.index[ fullTileNum, tileFragmentNum / 1k ] -> globalFragmentNumber (of tileFragmentNum-%1k)
 fragments.pos.index[ globalFragmentNumber / 1000 ] -> globalPos (of globalFragmentNumber-%1000)
 fragments.pos.shift[ globalFragmentNumber / 1000 ] -> position in file fragments.pos
 read&accumulate globalFragmentNumber%1000 entries from fragments.pos -> globalPos of globalFragmentNumber (of tileFragmentNum-%1k) { average 500 entries to go through }
 continue reading fragments.pos+fragments.tile, for tileFragmentNum%1k occurrences of fullTileNum -> globalPos of tileFragmentNum { average 500*256=128,000 entries to go through }

check v2:
 search in fragments.tile.index2 for closest smaller value to tileFragmentNum -> tileFragmentNum2 (<tileFragmentNum) and globalPos of tileFragmentNum2 (globalPos being a multiple of "a distance that'll contain 256,000 fragments" = 1.5Mbases at 30x) {11 comparisons, very fast}
 fragments.tile.index2[ "global", globalPos / 1000 ] -> globalFragmentNumber of tileFragmentNum2
 fragments.pos.shift2[ globalPos / 1000 ] -> position in file fragments.pos
 read fragments.pos+fragments.tile from there, for (tileFragmentNum-tileFragmentNum2) occurrences of fullTileNum -> globalPos of tileFragmentNum { average 128,000 entries to go through }


*/


        return 0;
    }


private:
};


} // namespace model
} // namespace eagle

#endif // EAGLE_MODEL_FRAGMENT_POS_RESOLVER_HH
