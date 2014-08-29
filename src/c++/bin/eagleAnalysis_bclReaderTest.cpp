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

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif //ifdef _OPENMP

using namespace std;


struct BclArguments
{
    string runFolder;
    string lane;
    string tile;
    string output;
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

void go( const BclArguments &args )
{
#ifdef _OPENMP
    omp_set_nested( 1 );
    //    omp_set_num_threads( 9 );
#endif //ifdef _OPENMP

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

        bclRunFolder.initReader( nextLane, nextTile, 1, cycleCount );

        while ( bclRunFolder.getRead( bclValues ) )
        {
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

            bclRunFolder.reportProgress( 0.1 );
        }
    }

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


    // Output
    clog << "Outputing main file..." << endl;
    ofstream outf( args.output.c_str(), ios_base::binary );
    outf.write( reinterpret_cast<char *>( &counts[0] ), counts.size() * sizeof( unsigned int ) );
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

int main( const int argc, const char **argv )
{
    return 0;
}

#endif //ifdef HAVE_LIBZOO
