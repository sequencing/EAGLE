/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "config.h"

#ifdef HAVE_LIBZOO

//#include "Logger.hh"
#include "libzoo/io/Bcl.hh"
#include "libzoo/cli/Common.hh"
#include "libzoo/util/AutoGrowVector.hh"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <boost/tuple/tuple.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif //ifdef _OPENMP

using namespace std;
using boost::tuple;
using boost::make_tuple;


AutoGrowVector< AutoGrowVector< AutoGrowVector< unsigned int > > > glitchStats; //[cycle][glitchHeight][glitchLength]);
AutoGrowVector< AutoGrowVector< AutoGrowVector< unsigned int > > > glitchStats2;
AutoGrowVector< unsigned int > glitchCountHist;
AutoGrowVector< unsigned int > glitchCountHist2;
vector< tuple< int, int, int > > glitchInfo2;
AutoGrowVector< AutoGrowVector< AutoGrowVector< unsigned int > > > qualityTable;


struct BclArguments
{
    string runFolder;
    string lane;
    string tile;
    string output;
    int firstCycle;
    int lastCycle;
    string verbosityLevel;
};

void printUsage()
{
    cout << "Parameters:" << endl;
    cout << "    --run-folder (-r)        Run folder path" << endl;
    cout << "    --lane (-l)              Lane filename (e.g. L001)" << endl;
    cout << "    --tile (-t)              Tile filename in lane (e.g s_1_1101)" << endl;
    cout << "    --output (-o)            Output prefix" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "    --verbose                = quiet [quiet|verbose|very-verbose|debug] or [0|1|2|3]" << endl;
    cout << "    --help (-h)              Help" << endl;
    cout << endl;
}

struct RepeatInfo
{
    unsigned int firstCycle;
    unsigned int lastCycle;
    unsigned int repeatKmerLength;
    unsigned int repeatCount;
    unsigned int repeatLength;

    RepeatInfo( const unsigned int i1, const unsigned int i2, const unsigned int i3, const unsigned int i4, const unsigned int i5 )
        : firstCycle( i1 )
        , lastCycle( i2 )
        , repeatKmerLength( i3 )
        , repeatCount( i4 )
        , repeatLength( i5 )
        {}
};

bool detectRepeatedKmer( const int maxKmerLength, const unsigned int repeatLengthThreshold, vector<unsigned int> &bclValues, unsigned int &startCycle, unsigned int &longestRepeat_endCycle, vector< RepeatInfo > &repeats )
{
    bool repeatDetected = false;
    unsigned int longestRepeat_repeatLength = 0;
    unsigned int longestRepeat_startCycle = 0;
    unsigned int longestRepeat_repeatCount = 0;
    int longestRepeat_repeatKmerLength = 0;

    for ( int repeatKmerLength = 1; repeatKmerLength <= maxKmerLength; ++repeatKmerLength )
    {
        unsigned int cycle = startCycle;

        while ( cycle + repeatKmerLength <= bclValues.size()
                && ( bclValues[cycle - 1] & 3 ) == ( bclValues[cycle - 1 + repeatKmerLength] & 3 ) )
        {
            ++cycle;
        }
        unsigned int endCycle = cycle + repeatKmerLength - 1;
        unsigned int repeatCount = ( endCycle - startCycle + 1 ) / repeatKmerLength;
        unsigned int repeatLength = ( repeatCount - 1 ) * repeatKmerLength; // Note: length excludes first element
        if ( repeatLength >= repeatLengthThreshold )
        {
            if ( repeatLength > longestRepeat_repeatLength )
            {
                repeatDetected = true;
                longestRepeat_repeatLength = repeatLength;
                longestRepeat_startCycle = startCycle;
                longestRepeat_endCycle = endCycle;
                longestRepeat_repeatCount = repeatCount;
                longestRepeat_repeatKmerLength = repeatKmerLength;
            }
        }
    }
    if ( repeatDetected )
    {
//        cout << "longestRepeat_Repeat found? startCycle=" << longestRepeat_startCycle << ", endCycle=" << longestRepeat_endCycle << ", repeatKmerLength=" << longestRepeat_repeatKmerLength << ", repeatCount=" << longestRepeat_repeatCount << ", repeatLength=" << longestRepeat_repeatLength << endl;
        repeats.push_back( RepeatInfo( longestRepeat_startCycle, longestRepeat_endCycle, longestRepeat_repeatKmerLength, longestRepeat_repeatCount, longestRepeat_repeatLength ) );
        startCycle = longestRepeat_endCycle;
    }
    return repeatDetected;
}

vector< RepeatInfo > detectFirstAreaWithRepeatsAndAllowedGap( const int maxKmerLength, const unsigned int repeatLengthThreshold, const unsigned int allowedGap, vector<unsigned int> &bclValues )
{
    vector< RepeatInfo > repeats;
    unsigned int longestRepeat_endCycle = 0;
    for ( unsigned int startCycle = 1;
          startCycle <= bclValues.size()
          && (
              longestRepeat_endCycle == 0 // If no repeat found, search until we find one
              || startCycle <= longestRepeat_endCycle + allowedGap // otherwise only attempt to continue the last found repeat with an allowed gap of 2-3 bases
          )
          ; ++startCycle )
    {
        detectRepeatedKmer( 10, 4, bclValues, startCycle, longestRepeat_endCycle, repeats );
    }
    return repeats;
}

void detectGlitchesInRead( vector<unsigned int> &bclValues, vector<int> &qualBinSteps, vector<int> &glitchInfo )
{
    unsigned int glitchCount = 0;
    for ( unsigned int cycle = 1; cycle <= bclValues.size(); ++cycle )
    {
        if ( qualBinSteps[cycle - 1] < -2 )
        {
            // maybe a new glitch
            int glitchHeight = -qualBinSteps[cycle - 1];
            int recoveredHeight = 0;
            int glitchLength = 0;
            // Glitch half life calculation
            while ( recoveredHeight < glitchHeight / 2
                    && cycle - 1 + glitchLength + 1 < qualBinSteps.size()
                    && qualBinSteps[cycle - 1 + glitchLength + 1] >= 0 )
            {
                glitchLength++;
                recoveredHeight += qualBinSteps[cycle - 1 + glitchLength];
            }
            int glitchHalfLife = glitchLength;
            // Glitch almost full recovery
            while ( recoveredHeight < glitchHeight - 1
                    && cycle - 1 + glitchLength + 1 < qualBinSteps.size()
                    && qualBinSteps[cycle - 1 + glitchLength + 1] >= 0 )
            {
                glitchLength++;
                recoveredHeight += qualBinSteps[cycle - 1 + glitchLength];
            }
            // Glitch full recovery
            if ( recoveredHeight < glitchHeight
                 && cycle - 1 + glitchLength + 1 < qualBinSteps.size()
                 && qualBinSteps[cycle - 1 + glitchLength + 1] > 0 )
            {
                glitchLength++;
                recoveredHeight += qualBinSteps[cycle - 1 + glitchLength];
            }

            if ( recoveredHeight >= glitchHeight / 2
                 && glitchLength < 8
               )
            {
                ++glitchCount;
                glitchInfo[cycle - 1] = -1;
                int i;
                for ( i = 1; i <= glitchHalfLife; i++ )
                    glitchInfo[cycle - 1 + i] = 1;
                for ( ; i <= glitchLength; i++ )
                    glitchInfo[cycle - 1 + i] = 2;

                // Add glitch to stats
                ++(glitchStats[cycle][glitchHeight][glitchLength]);
                glitchInfo2.push_back(make_tuple(cycle,glitchHeight,glitchLength));
            }
        }
    }
    ++(glitchStats[0][0][0]);
    ++(glitchCountHist[glitchCount]);
}

void go( const BclArguments &args )
{
    // Check that we won't overwrite any existing file
    if ( doesFileExist( args.output ) )
    {
        cerr << args.output << " already exists in the current directory. Aborting." << endl;
        exit( -1 );
    }


    BclRunFolder bclRunFolder( args.runFolder, args.lane, args.tile );
    unsigned int cycleCount = bclRunFolder.getCycleCount();
    //    cycleCount /= 2;
    //    --cycleCount;
    //    clog << "Processing read 1 only: " << cycleCount << " cycles" << endl;
    int firstCycle = args.firstCycle;
    int lastCycle = args.lastCycle;
    if (firstCycle == 0)
        firstCycle = 1;
    if (lastCycle == 0)
        lastCycle = cycleCount / 2;

    cycleCount = lastCycle - firstCycle + 1;

    vector<unsigned int> bclValues;
    vector<unsigned int> counts( cycleCount * 1024 * 4 * 8 * 8 );
    int Qbins[41] =
    {
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
        3, 3, 3, 3, 3, 4, 4, 4, 4, 4,
        4, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7
    };


    string nextLane, nextTile;
    while ( bclRunFolder.getNextLaneAndTileNames( nextLane, nextTile ) )
    {
        bclRunFolder.initReader( nextLane, nextTile, firstCycle, lastCycle );

        while ( bclRunFolder.getRead( bclValues ) )
        {
            vector<int> qualBins;
            vector<int> qualBinSteps;
            //qualBinSteps.clear();

            for ( unsigned int cycle = 1; cycle <= bclValues.size(); ++cycle )
            {
                qualBins.push_back( Qbins[bclValues[cycle - 1] >> 2] );
            }
            for ( unsigned int cycle = 1; cycle <= bclValues.size() - 1; ++cycle )
            {
                qualBinSteps.push_back( qualBins[cycle] - qualBins[cycle - 1] );
            }

            // Glitch detection
            vector<int> glitchInfo( bclValues.size(), 0 );
            glitchInfo2.clear();
            detectGlitchesInRead( bclValues, qualBinSteps, glitchInfo );

            // Motif detection: Detect repetitive 1to10mers, with minimum size 4 and allowed gap 3
            vector< RepeatInfo > repeats = detectFirstAreaWithRepeatsAndAllowedGap( 10, 4, 20, bclValues );
            if ( repeats.size() == 1)
            {
                unsigned int first = repeats.front().firstCycle;
                unsigned int last = repeats.back().lastCycle;

                // Calculate average Q in the next 20 cycles
                unsigned int averageQ = 0;
                unsigned int countAfter = min<unsigned int>( bclValues.size() - last, 20 );
                unsigned int countAfterExcludingGlitches = 0;
                for ( unsigned int cycle = last+1; cycle <= last+countAfter; ++cycle )
                {
                    bool excludeGlitches = false;
                    if (!excludeGlitches || glitchInfo[cycle - 1] == 0)
                    {
                        averageQ += bclValues[cycle - 1] >> 2;
                        ++countAfterExcludingGlitches;
                    }
                }
                if (countAfterExcludingGlitches > 5)
                {
                    averageQ /= countAfterExcludingGlitches;

                    // Find smallest kmer rotation
                    {
                        unsigned long long kmer = 0;
                        for ( unsigned int i = first; i < first+repeats.front().repeatKmerLength; ++i )
                        {
                            kmer <<= 2;
                            kmer |= bclValues[i - 1] & 3;
                        }
                        unsigned long long kmerMask = (1 << (2*repeats.front().repeatKmerLength)) - 1;
                        unsigned long long smallestRotatedKmer = kmer;
//                        unsigned int smallestRot = 0;
                        for (unsigned int rot = 1; rot < repeats.front().repeatKmerLength; ++rot )
                        {
                            unsigned long long rotatedBase = kmer >> (2*repeats.front().repeatKmerLength - 2);
                            kmer <<= 2;
                            kmer &= kmerMask;
                            kmer |= rotatedBase;
                            if (kmer < smallestRotatedKmer)
                            {
                                smallestRotatedKmer = kmer;
//                                smallestRot = rot;
                            }
                        }

                        kmer = smallestRotatedKmer;
                        for ( unsigned int i = 0; i < repeats.front().repeatKmerLength; ++i )
                        {
                            unsigned long long lhsBase = kmer >> (2*repeats.front().repeatKmerLength - 2);
                            kmer <<= 2;
                            kmer &= kmerMask;
                            cout << lhsBase;
                        }
                    }

                    cout << "\t" << repeats.front().repeatCount;
                    cout << "\t" << ( last - first + 1 ); //( repeats.front().repeatLength + repeats.front().repeatKmerLength );
                    cout << "\t" << averageQ;
                    cout << "\t" << first << "->" << last;

                    cout << "\t";
                    for ( unsigned int i = first; i < first+repeats.front().repeatKmerLength; ++i )
                    {
                        cout << ( bclValues[i - 1] & 3 );
                    }

                    cout << "\t" << repeats.front().repeatKmerLength;

                    cout << "\t";
                    for ( unsigned int i = first; i <= last; ++i )
                    {
                        cout << ( bclValues[i - 1] & 3 );
                    }

                    cout << "\t";
                    for ( unsigned int i = 1; i <= bclValues.size(); ++i )
                    {
                        if ( i == first )
                            cout << " | ";
                        cout << ( bclValues[i - 1] & 3 );
                        if ( i == last )
                            cout << " | ";
                    }
                    cout << endl;
                }
            }
            else if ( repeats.size() == 0)
            {
                // No repeat found => output glitch to the safe list
                for (unsigned int i=0; i<glitchInfo2.size(); ++i)
                    ++(glitchStats2
                       [glitchInfo2[i].get<0>()]
                       [glitchInfo2[i].get<1>()]
                       [glitchInfo2[i].get<2>()]
                        );
                ++(glitchStats2[0][0][0]);
                ++(glitchCountHist2[glitchInfo2.size()]);
            }

            // Final quality table stats
            if ( repeats.size() == 0) // Let's try to use only those perfect reads at first
            {
                int averageQ = 0;
                int qCount = 0;
                for ( unsigned int cycle = 1; cycle <= bclValues.size(); ++cycle )
                {
                    if (glitchInfo[cycle - 1] == 0)
                    {
                        averageQ += bclValues[cycle - 1] >> 2;
                        ++qCount;
                    }
                }
                if (qCount > 20)
                {
                    averageQ /= qCount;
                    int profileId = averageQ;
                    ++(qualityTable[0][0][profileId]);
                    for ( unsigned int cycle = 1; cycle <= bclValues.size(); ++cycle )
                    {
                        if (glitchInfo[cycle - 1] == 0)
                        {
                            int qual = bclValues[cycle - 1] >> 2;
                            ++(qualityTable[profileId][cycle][qual]);
                        }
                    }
                }
            }


            // old stuff
/*
            unsigned int precedingKmer = 0;
            unsigned int prevQualityBin = 0;
            for ( unsigned int cycleNum = 1; cycleNum <= cycleCount; ++cycleNum )
            {
                if ( cycleNum > 1 )
                {
                    precedingKmer = ( ( precedingKmer & 0xFF ) << 2 ) | ( bclValues[cycleNum - 2] & 3 );
                    prevQualityBin = Qbins[bclValues[cycleNum - 2] >> 2];
                }
                unsigned int newBaseNum = bclValues[cycleNum - 1] & 3;
                unsigned int newQualityBin = Qbins[bclValues[cycleNum - 1] >> 2];

                assert( cycleNum >= 1 && cycleNum <= cycleCount );
                assert( precedingKmer < 1024 );
                assert( newBaseNum < 4 );
                assert( prevQualityBin < 8 );
                assert( newQualityBin < 8 );
                unsigned int entryNum =
                    ( cycleNum - 1 ) * ( 1024 * 4 * 8 * 8 ) +
                    precedingKmer * ( 4 * 8 * 8 ) +
                    newBaseNum * ( 8 * 8 ) +
                    prevQualityBin * ( 8 ) +
                    newQualityBin;

                ++( counts[entryNum] ) ;
            }
*/

            bclRunFolder.reportProgress( 0.1 );
        }
    }

/*
    // Report
    clog << "Calculating report..." << endl;
    unsigned int sum = 0;
    unsigned int count = 0;
    unsigned int maxValue = 0;
    unsigned int maxIndex = 0;
    for ( unsigned i = 0;  i < counts.size(); ++i )
    {
        if ( counts[i] )
        {
            ++count;
            sum += counts[i];
            if ( ( i > 5 * 1024 * 4 * 8 * 8 ) && ( maxValue < counts[i] ) )
            {
                maxValue = counts[i];
                maxIndex = i;
            }
        }
    }
    clog << "count=" << count << " = " << count * 100.0 / counts.size() << "%" << endl;
    clog << "sum=" << sum << endl;
    clog << "max=" << maxValue << " for index=" << maxIndex
         << "={ cycle=" << maxIndex / ( 1024 * 4 * 8 * 8 )
         << ", precedingKmer=" << maxIndex / ( 4 * 8 * 8 ) % 1024
         << ", newBase=" << maxIndex / ( 8 * 8 ) % 4
         << ", prevQbin=" << maxIndex / ( 8 ) % 8
         << ", newQbin=" << maxIndex / ( 1 ) % 8
         << " }" << endl;
*/


    // Output
    clog << "Outputing glitch file..." << endl;
    ofstream outf( (args.output + ".glitches").c_str(), ios_base::binary );
//    outf.write( reinterpret_cast<char *>( &counts[0] ), counts.size() * sizeof( unsigned int ) );
    for (unsigned int i=0; i<glitchStats.size(); ++i)
        for (unsigned int j=0; j<glitchStats[i].size(); ++j)
            for (unsigned int k=0; k<glitchStats[i][j].size(); ++k)
                outf << i << "\t" << j << "\t" << k << "\t" << glitchStats[i][j][k] << endl;

    ofstream outf2( (args.output + ".glitches2").c_str(), ios_base::binary );
    for (unsigned int i=0; i<glitchStats2.size(); ++i)
        for (unsigned int j=0; j<glitchStats2[i].size(); ++j)
            for (unsigned int k=0; k<glitchStats2[i][j].size(); ++k)
                outf2 << i << "\t" << j << "\t" << k << "\t" << glitchStats2[i][j][k] << endl;

    ofstream outf3( (args.output + ".glitchCountHist").c_str(), ios_base::binary );
    for (unsigned int i=0; i<glitchCountHist.size(); ++i)
        outf3 << i << "\t" << glitchCountHist[i] << endl;

    ofstream outf4( (args.output + ".glitchCountHist2").c_str(), ios_base::binary );
    for (unsigned int i=0; i<glitchCountHist2.size(); ++i)
        outf4 << i << "\t" << glitchCountHist2[i] << endl;


    clog << "Outputing quality table..." << endl;
    ofstream outQTable( (args.output + ".qtable").c_str(), ios_base::binary );
//    outf.write( reinterpret_cast<char *>( &counts[0] ), counts.size() * sizeof( unsigned int ) );
    for (unsigned int i=0; i<qualityTable.size(); ++i)
        for (unsigned int j=0; j<qualityTable[i].size(); ++j)
            for (unsigned int k=0; k<qualityTable[i][j].size(); ++k)
                if (qualityTable[i][j][k])
                    outQTable << i << "\t" << j << "\t" << k << "\t" << qualityTable[i][j][k] << endl;

    ofstream outQTable2( (args.output + ".qtable2").c_str(), ios_base::binary );
//    outf.write( reinterpret_cast<char *>( &counts[0] ), counts.size() * sizeof( unsigned int ) );
    for (unsigned int i=0; i<qualityTable.size(); ++i)
        for (unsigned int j=0; j<qualityTable[i].size(); ++j)
        {
            bool dataFound = false;
            for (unsigned int k=0; k<qualityTable[i][j].size(); ++k)
                if (qualityTable[i][j][k])
                {
                    dataFound = true;
                    break;
                }

            if (dataFound)
            {
                outQTable2 << i << "\t" << j;
                for (unsigned int k=0; k<qualityTable[i][j].size(); ++k)
                    if (qualityTable[i][j][k])
                        outQTable2 << "\t" << k << ":" << qualityTable[i][j][k];
                outQTable2 << endl;
            }
        }
}
/*
    {
        cout << "Read:";
        for (unsigned int i=0; i<bclValues.size(); ++i)
        {
            cout << (int)bclValues[i] << ",";
        }
        cout << endl;
*/


int main( const int argc, const char **argv )
{
    cout << "Command called:" << endl << "   ";
    for ( int i = 0; i < argc; ++i )
    {
        cout << " " << argv[i];
    }
    cout << "\n" << endl;

    BclArguments args;

    for ( int i = 1; i < argc; ++i )
    {
        if ( isNextArgument( "-h", "--help", argc, argv, i, NULL ) )
        {
            printUsage();
            exit( 0 );
        }
        else if ( isNextArgument   ( "-r", "--run-folder"          , argc, argv, i, &args.runFolder                       ) ) {}
        else if ( isNextArgument   ( "-l", "--lane"                , argc, argv, i, &args.lane                            ) ) {}
        else if ( isNextArgument   ( "-t", "--tile"                , argc, argv, i, &args.tile                            ) ) {}
        else if ( isNextArgument   ( "-o", "--output"              , argc, argv, i, &args.output                          ) ) {}
        else if ( isNextArgumentInt( "-b", "--first-cycle"         , argc, argv, i, &args.firstCycle                      ) ) {}
        else if ( isNextArgumentInt( "-e", "--last-cycle"          , argc, argv, i, &args.lastCycle                       ) ) {}
        else if ( isNextArgument   ( ""  , "--verbose"             , argc, argv, i, &args.verbosityLevel                  ) )
        {
            //            Logger::setVerbosity( args.verbosityLevel );
        }
        else
        {
            cerr << "Error: Invalid parameter: " << argv[i] << "\n" << endl;
            printUsage();
            exit( 1 );
        }
    }

    // Checking for required parameters
    if ( args.runFolder.empty() || args.output.empty() )
    {
        cerr << "Error: Missing arguments: --run-folder, --lane, --tile and --output are required\n" << endl;
        printUsage();
        exit( 1 );
    }

    // Go
    go( args );

    return 0;
}


#else //ifdef HAVE_LIBZOO

#include <iostream>

int main( const int argc, const char **argv )
{
    std::cerr << "EAGLE must be compiled with libzoo to get this tool" << std::endl;
    return 0;
}

#endif //ifdef HAVE_LIBZOO
