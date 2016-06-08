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

#include "libzoo/cli/Common.hh"
#include "libzoo/io/Bcl.hh"

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _OPENMP
#include <omp.h>
#endif //ifdef _OPENMP

using namespace std;


std::vector< unsigned int > bigTable_;

void importBigQualityTableFile( const string &filename )
{
    struct stat fileStat;
    stat( filename.c_str(), &fileStat );
    unsigned long totalSize = fileStat.st_size;
    bigTable_.resize( totalSize / sizeof( unsigned int ) );
    ifstream is( filename.c_str(), ios_base::binary );
    is.read( reinterpret_cast<char *>( &bigTable_[0] ), totalSize );
}

/* not used yet
static int Qbins[41] =
{
    0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
    3, 3, 3, 3, 3, 4, 4, 4, 4, 4,
    4, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7
};

static int bin2Q[8] = { 0, 6, 15, 22, 27, 33, 37, 40 };
*/

void printKmer( unsigned int kmer, int length )
{
    for ( int i = 0; i < length; ++i )
    {
        switch ( ( kmer >> ( 2 * length - 2 ) ) & 3 )
        {
            case 0:
                cout << 'A';
                break;
            case 1:
                cout << 'C';
                break;
            case 2:
                cout << 'G';
                break;
            case 3:
                cout << 'T';
                break;
        }
        kmer <<= 2;
    }
}


unsigned int getQuality( const int cycle, int precedingKmer, int newBaseNum, const int prevQBin )
{
    //    unsigned int prevQualityBin = Qbins[ clusterErrorModelContext.qualityModelContext.profileNumber ];
    unsigned int firstCycle = ( cycle > 0 ) ? cycle : 1;
    unsigned int lastCycle = ( cycle > 0 ) ? cycle : 502;
    unsigned int firstPrevQBin = ( prevQBin >= 0 ) ? prevQBin : 0;
    unsigned int lastPrevQBin = ( prevQBin >= 0 ) ? prevQBin : 7;
    int firstPrecedingKmer = ( precedingKmer >= 0 ) ? precedingKmer : 0;
    int lastPrecedingKmer = ( precedingKmer >= 0 ) ? precedingKmer : 0x3FF;
    int firstNewBaseNum = ( newBaseNum >= 0 ) ? newBaseNum : 0;
    int lastNewBaseNum = ( newBaseNum >= 0 ) ? newBaseNum : 3;
    for ( precedingKmer = firstPrecedingKmer; precedingKmer <= lastPrecedingKmer; ++precedingKmer )
    {
        for ( newBaseNum = firstNewBaseNum; newBaseNum <= lastNewBaseNum; ++newBaseNum )
        {
            vector<unsigned int> counts( 8, 0 );
            unsigned int sum = 0;
            for ( unsigned int cycleNum = firstCycle; cycleNum <= lastCycle; ++cycleNum )
            {
                for ( unsigned int prevQualityBin = firstPrevQBin; prevQualityBin <= lastPrevQBin; ++prevQualityBin )
                {
                    for ( unsigned int newQualityBin = 0; newQualityBin < 8; ++newQualityBin )
                    {
                        assert( cycleNum >= 1 && cycleNum <= 502 );
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
                        assert( entryNum < bigTable_.size() );
                        unsigned int val = bigTable_[entryNum];
                        counts[newQualityBin] += val;
                        sum += val;
                    }
                }
            }
            cout << "sum=" << sum << endl;

            if ( sum != 0 )
            {
                for ( unsigned int newQualityBin = 0; newQualityBin < 8; ++newQualityBin )
                {
                    cout << "kmer=";
                    printKmer( precedingKmer, 5 );
                    cout << "\tnewBase=";
                    printKmer( newBaseNum, 1 );
                    cout << "\tnewQ=" << newQualityBin << "\tcount=" << counts[newQualityBin] << "\t" << counts[newQualityBin] * 100.0 / sum << endl;
                }
            }
        }
    }

    return 0;
}


struct BclArguments
{
    string table;
    string kmer;
    string cycle;
    string prevQBin;
    string newBase;
    string verbosityLevel;
};

void printUsage()
{
    cout << "Parameters:" << endl;
    cout << "    --table (-t)             " << endl;
    cout << "    --kmer (-k)              " << endl;
    cout << "    --cycle (-c)              " << endl;
    cout << "    --prev-qbin (-p)              " << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "    --verbose                = quiet [quiet|verbose|very-verbose|debug] or [0|1|2|3]" << endl;
    cout << "    --help (-h)              Help" << endl;
    cout << endl;
}

unsigned int kmerAsciiToNum( string s )
{
    if ( ( s[0] >= '0' && s[0] <= '9' ) || s[0] ==  '-' )
        return atoi( s.c_str() );

    unsigned int kmer = 0;
    for ( unsigned i = 0; i < s.size(); ++i )
    {
        kmer <<= 2;
        switch ( s[i] )
        {
            case 'A':
                kmer += 0;
                break;
            case 'C':
                kmer += 1;
                break;
            case 'G':
                kmer += 2;
                break;
            case 'T':
                kmer += 3;
                break;
            default:
                assert( false );
        }
    }
    return kmer;
}

void go( const BclArguments &args )
{
    importBigQualityTableFile( args.table );
    const unsigned int precedingKmer = kmerAsciiToNum( args.kmer );
    const unsigned int cycleNum = atoi( args.cycle.c_str() );
    const unsigned int prevQualityBin = atoi( args.prevQBin.c_str() );
    const unsigned int newBaseNum = kmerAsciiToNum( args.newBase );
    getQuality( cycleNum, precedingKmer, newBaseNum, prevQualityBin );
}


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
        else if ( isNextArgument   ( "-t", "--table"          , argc, argv, i, &args.table                       ) ) {}
        else if ( isNextArgument   ( "-k", "--kmer"                , argc, argv, i, &args.kmer                            ) ) {}
        else if ( isNextArgument   ( "-c", "--cycle"                , argc, argv, i, &args.cycle                            ) ) {}
        else if ( isNextArgument   ( "-p", "--prev-qbin"                , argc, argv, i, &args.prevQBin                            ) ) {}
        else if ( isNextArgument   ( "-b", "--new-base"                , argc, argv, i, &args.newBase                            ) ) {}
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
    /*
        if ( args.runFolder.empty() || args.output.empty() )
        {
            cerr << "Error: Missing arguments: --run-folder, --lane, --tile and --output are required\n" << endl;
            printUsage();
            exit( 1 );
        }
    */

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
