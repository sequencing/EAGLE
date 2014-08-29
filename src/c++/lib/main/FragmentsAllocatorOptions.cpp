/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Command line options for 'allocateFragments'
 **
 ** \author Lilian Janin
 **/

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "main/FragmentsAllocatorOptions.hh"

namespace eagle
{
namespace main
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

FragmentsAllocatorOptions::FragmentsAllocatorOptions()
    : sampleGenomeDir( bfs::current_path() / "sample_genome" )
    , outputDir( bfs::current_path() )
    , coverageDepth( 30 )
    , tls( "380:400:420:10:10:FRp:RFm" )
    , uniformCoverage( false )
    , tileAllocationMethodStr( "random" )
    , randomSeed( 1 )
    , templateLengthTableFile("")
    , contigName( "" )
    , mergeExistingFragments( false )
    , gcCoverageFitFile( "" )
{
    parameters_.add_options()
        ("sample-genome-dir,s",      bpo::value< bfs::path >(&sampleGenomeDir)->default_value( sampleGenomeDir, "./sample_genome" ),
                                     "[input]  \tFull path to the directory containing the sample's genome FASTA files")
        ("output-dir,f",             bpo::value< bfs::path >(&outputDir)->default_value( outputDir, "." ),
                                     "[output] \tFull path to the location where the fragments should be written")
        ;
    namedOptions_.add_options()
        ("coverage-depth,d",         bpo::value< float >(&coverageDepth)->default_value(coverageDepth),
                                     "Desired coverage depth")
        ("tiles,t",                  bpo::value< unsigned long >(&tileCount),
                                     "Number of desired tiles")
        ("bases-per-cluster,b",      bpo::value< unsigned int >(&basesPerCluster),
                                     "Number of bases per cluster of reads (i.e. cluster length)")
        //        ("reads,r", bpo::value< unsigned long >(&readCount)->default_value(1000000), "Number of desired reads")
        ("tls",                      bpo::value< std::string >(&tls)->default_value( tls ),
                                     "Template-length statistics in the format 'min:median:max:lowStdDev:highStdDev:M0:M1', "
                                     "where M0 and M1 are the numeric value of the models (0=FFp, 1=FRp, 2=RFp, 3=RRp, 4=FFm, 5=FRm, 6=RFm, 7=RRm) "
                                     "- only min and max are currently used, but all are kept to match iSAAC format")
        ("uniform-coverage,u",       bpo::value< bool >(&uniformCoverage)->zero_tokens(),
                                     "Generates equally-spaced reads across all tiles with a fixed template length equal to the specified median")
        ("tile-allocation-method",   bpo::value< std::string >(&tileAllocationMethodStr)->default_value(tileAllocationMethodStr),
                                     "Possible values are: random, sequence (we start filling tile 1 with reads until full, then tile 2, etc.), "
                                     "interleaved (read 1 goes to tile 1, read 2 -> tile 2, ..., read N -> tile N, read N+1 -> tile 1, etc.)")
        ("random-seed",              bpo::value<unsigned int>(&randomSeed)->default_value(randomSeed),
                                     "Seed to use for random number generation")
        ("template-length-table",    bpo::value< bfs::path >(&templateLengthTableFile),
                                     "File containing the template length table")
        ("contig",                   bpo::value< std::string >(&contigName)->default_value(contigName),
                                     "If specified, only generate fragments for this contig")
        ("merge-existing-fragments", bpo::value< bool >(&mergeExistingFragments)->zero_tokens(),
                                     "Merge pre-calculated fragment files. Don't compute new fragments.")
        ("gc-coverage-fit-table",    bpo::value< bfs::path >(&gcCoverageFitFile),
                                     "File describing how GC content affects the coverage")
        ("max-coverage-error",       bpo::value< double >(&maxCoverageError)->default_value( 0.25 ),
                                     "Using --gc-coverage-fit can lead to inexact coverage depth. If the error is over this threshold, it will restart and try to do better")
        ;
}

void FragmentsAllocatorOptions::postProcess(bpo::variables_map &vm)
{
    eagle::common::OptionsHelper check(vm);

    check.requiredOptions(boost::assign::list_of
                          ("coverage-depth")
                          ("bases-per-cluster")
                          ("tiles")
                          );

    // Parsing "Template-length statistics"
    std::vector<std::string> tokens;
    boost::split( tokens, tls, boost::is_any_of(":"));
    if ( tokens.size() == 7 )
    {
        templateLengthStatistics.min        = boost::lexical_cast<double>( tokens[0] );
        templateLengthStatistics.median     = boost::lexical_cast<double>( tokens[1] );
        templateLengthStatistics.max        = boost::lexical_cast<double>( tokens[2] );
        templateLengthStatistics.lowStdDev  = boost::lexical_cast<double>( tokens[3] );
        templateLengthStatistics.highStdDev = boost::lexical_cast<double>( tokens[4] );
        templateLengthStatistics.M0         = tokens[5];
        templateLengthStatistics.M1         = tokens[6];
    }
    else
    {
        EAGLE_ERROR( "--tls option must have 7 colon-separated tokens" );
    }

    // Parsing "Tile allocation method"
    if      (tileAllocationMethodStr == "random"     ) { tileAllocationMethod = TILE_ALLOCATION_RANDOM; }
    else if (tileAllocationMethodStr == "sequence"   ) { tileAllocationMethod = TILE_ALLOCATION_SEQUENCE; }
    else if (tileAllocationMethodStr == "interleaved" ||
             tileAllocationMethodStr == "interleave" ) { tileAllocationMethod = TILE_ALLOCATION_INTERLEAVED; }
    else
    {
        EAGLE_ERROR( "Invalid value for --tile-allocation-method option" );
    }
}

} //namespace main
} // namespace eagle
