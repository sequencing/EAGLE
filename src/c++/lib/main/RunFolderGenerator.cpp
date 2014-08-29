/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include <boost/assign.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/integer.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <fstream>
#include <cmath>

#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "main/RunFolderGenerator.hh"

using namespace std;


namespace eagle
{
namespace main
{

    /*
    const unsigned int LANE_COUNT = 8;

    const unsigned int READ_LENGTH = 101;
    const unsigned int BARCODE_LENGTH = 7;
    const unsigned int CLUSTER_LENGTH = 2 * READ_LENGTH;

    const std::vector<unsigned int> lengths = boost::assign::list_of(READ_LENGTH)(BARCODE_LENGTH)(READ_LENGTH);
    const std::string instrument = "EAS987";
    */
    /*
    using boost::posix_time::to_iso_string;
    const std::string today = to_iso_string(boost::gregorian::day_clock::local_day());
    const std::string date = today.substr(2);
    */
    /*
    const std::string flowcell = "FC1234TST";
    const unsigned int runFolderId = 567;
    */

/* We now provide tile names as a command line parameter
    const std::vector<unsigned int> tileNameMethod1 = boost::assign::list_of(1)(2)(3)(4)(5)(6)(7)(8)
        (21)(22)(23)(24)(25)(26)(27)(28)
        (41)(42)(43)(44)(45)(46)(47)(48)
        (61)(62)(63)(64)(65)(66)(67)(68);
*/

    //    const std::vector<std::string> barcodes;// = boost::assign::list_of("ACACACA");




RunFolderGenerator::RunFolderGenerator( const RunFolderGeneratorOptions &options )
    : options_        ( options )
    , runInfo_        ( options_.runInfo )
    , runFolderPath_  ( options_.outputDir )
    , dataPath_       ( runFolderPath_ / "Data" )
    , intensitiesPath_( dataPath_ / "Intensities" )
    , baseCallsPath_  ( intensitiesPath_ / "BaseCalls" )
{
    std::clog << "+ outputDir: " << options_.outputDir << std::endl;

    runFolder_ = runInfo_.runId;
    //    const std::string runFolder = (boost::format("%s_%s_%04d_%s") % date % instrument % runId % flowcell).str();
    vector<string> tok;
    boost::split( tok, runFolder_, boost::is_any_of("_") );
    assert( tok.size() == 4 );
    runFolderDate_ = tok[0];
    instrument_ = tok[1];
    runFolderId_ = tok[2];
    flowcell_ = tok[3];
}


void RunFolderGenerator::run()
{
    generateDirectoryStructure();
    generateMetadata();
}


void RunFolderGenerator::generateDirectoryStructure() const
{
    std::clog << "Generating directory structure " << runFolderPath_.string() << std::endl;

    create_directory( runFolderPath_ );
    create_directory( dataPath_ );
    create_directory( intensitiesPath_ );
    create_directory( baseCallsPath_ );

    for (unsigned int i=1; i<=runInfo_.laneCount; ++i)
    {
        std::clog << "    ... for lane " << i << std::endl;
        const bfs::path lanePath( baseCallsPath_ / (boost::format("L%03g") % i).str() );
        create_directory( lanePath );

        for (unsigned int j=1; j<=runInfo_.getClusterLength(); ++j)
        {
            const bfs::path clusterPath( lanePath / (boost::format("C%d.1") % j).str() );
            create_directory( clusterPath );
        }

        // Directory for locs/clocs files
        const bfs::path lanePath2( intensitiesPath_ / (boost::format("L%03g") % i).str() );
        create_directory( lanePath2 );
    }
}


void RunFolderGenerator::generateMetadata() const
{
    generateRunInfo();
    generateConfig();
    //    generateSampleSheet();
    generateMatrix();
    generatePhasing();
}

/*****************************************************************************/

void RunFolderGenerator::generateRunInfo() const
{
    using boost::property_tree::ptree;
    ptree pt;
    pt.put("RunInfo.<xmlattr>.xmlns:xsd", "http://www.w3.org/2001/XMLSchema");
    pt.put("RunInfo.<xmlattr>.xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    pt.put("RunInfo.<xmlattr>.Version", "2");
    pt.put("RunInfo.Run.<xmlattr>.Id", runFolder_);
    pt.put("RunInfo.Run.<xmlattr>.Number", runInfo_.runNumber);
    pt.put("RunInfo.Run.<xmlattr>.TileNameMethod", runInfo_.tileNameMethod);
    ptree &run = pt.get_child("RunInfo.Run");
    run.put("Flowcell", flowcell_);
    run.put("Date", runFolderDate_);
    run.put("FlowcellLayout.<xmlattr>.LaneCount", runInfo_.laneCount);
    run.put("FlowcellLayout.<xmlattr>.SurfaceCount", runInfo_.surfaceCount);
    run.put("FlowcellLayout.<xmlattr>.SwathCount", runInfo_.swathCount);
    run.put("FlowcellLayout.<xmlattr>.TileCount", runInfo_.tileCount);
    run.put("Reads", "");
    ptree &reads = run.get_child("Reads");

    int readNum = 1;
    BOOST_FOREACH(const eagle::io::ReadDescription &rd, runInfo_.reads)
    {
        ptree read;

        // Attributes for recent RunInfo.xml format
        read.put("<xmlattr>.Number", readNum++);
        read.put("<xmlattr>.IsIndexedRead", rd.isIndex?"Y":"N");
        read.put("<xmlattr>.NumCycles", rd.lastCycle - rd.firstCycle + 1);

        // Attributes for older RunInfo.xml format
        read.put("<xmlattr>.FirstCycle", rd.firstCycle);
        read.put("<xmlattr>.LastCycle", rd.lastCycle);
        if (rd.isIndex)
        {
            read.put("Index", "");
        }

        reads.add_child("Read", read);
    }

    const boost::property_tree::xml_writer_settings<char> w(' ', 2);
    const bfs::path runInfoPath = runFolderPath_ / "RunInfo.xml";

    // DEBUG
    std::clog << "Writing RunInfo file " << runInfoPath.string() << std::endl;

    write_xml(runInfoPath.string(), pt, std::locale(), w);
}

/*****************************************************************************/

void RunFolderGenerator::generateConfig() const
{
    using boost::property_tree::ptree;
    ptree pt;
    pt.put("BaseCallAnalysis.Run.<xmlattr>.Name", "BaseCalls");
    ptree &run = pt.get_child("BaseCallAnalysis.Run");
    run.put("Cycles.<xmlattr>.First", 1);
    const unsigned int length = runInfo_.getClusterLength();
    run.put("Cycles.<xmlattr>.Last", length);
    run.put("Cycles.<xmlattr>.Length", length);
    run.put("BaseCallParameters", "");
    ptree &baseCallParameters = run.get_child("BaseCallParameters");
    run.put("RunParameters", "");
    ptree &runParameters = run.get_child("RunParameters");

    baseCallParameters.put("ChastityThreshold", 0.6); // This entry is only needed to keep the structure of the xml file such that copyConfig.pl doesn't "optimise" things and generate a broken output file

    unsigned int r = 0;
    BOOST_FOREACH(const eagle::io::ReadDescription &rd, runInfo_.reads)
    {
        r++;
        ptree matrix;
        matrix.put("Read", r);
        matrix.put("AutoFlag", 2);
        matrix.put("AutoLane", 0);
        matrix.put("FirstCycle", rd.firstCycle);
        matrix.put("LastCycle", rd.lastCycle);
        ptree phasing = matrix;
        matrix.put("CycleOffset", 0);
        matrix.put("Cycle", rd.firstCycle);
        phasing.put("CycleOffset", 1);
        phasing.put("Cycle", rd.firstCycle + 1);
        phasing.put("PhasingRate", 0.02 + r / 100.0);
        phasing.put("PrephasingRate", 0.03 + r / 100.0);
        baseCallParameters.add_child("Matrix", matrix);
        baseCallParameters.add_child("Phasing", phasing);
        ptree reads;
        reads.put("<xmlattr>.Index", r);
        reads.put("FirstCycle", rd.firstCycle);
        reads.put("LastCycle", rd.lastCycle);
        runParameters.add_child("Reads", reads);
    }

    runParameters.put("Instrument", instrument_);
    runParameters.put("RunFolder", runFolder_);
    runParameters.put("RunFolderDate", runFolderDate_);
    runParameters.put("RunFolderId", runFolderId_);
    run.put("Software", "");
    run.put("Software.<xmlattr>.Name", "RTA");
    run.put("Software.<xmlattr>.Version", "1.9.35.0");
    run.put("TileSelection", "");
    ptree &tileSelection = run.get_child("TileSelection");

    for (unsigned int l = 1; runInfo_.laneCount >= l; ++l)
    {
        ptree lane;
        lane.put("<xmlattr>.Index", l);
        lane.put("Sample", 's');
        unsigned int tilesPerLane = runInfo_.tileCount * runInfo_.surfaceCount * runInfo_.swathCount;
        assert( tilesPerLane == options_.tileId.size() );
        for (unsigned int t = 0; t < tilesPerLane; ++t)
        {
            lane.add("Tile", options_.tileId[t]);
        }
        tileSelection.add_child("Lane", lane);
    }

    const boost::property_tree::xml_writer_settings<char> w(' ', 2);
    const bfs::path configPath = baseCallsPath_ / "config.xml";

    // DEBUG
    std::clog << "Writing config file " << configPath.string() << std::endl;

    write_xml(configPath.string(), pt, std::locale(), w);
}

/*****************************************************************************/
/*
void RunFolderGenerator::generateSampleSheet()
{
    const bfs::path sampleSheet = baseCallsPath_ / "SampleSheet.csv";

    // DEBUG
    std::clog << "Writing sample sheet " << sampleSheet.string() << std::endl;

    std::ofstream os(sampleSheet.string().c_str());
    os << "#FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject" << std::endl;

    for (unsigned int lane = 1; 8 >= lane; ++lane)
    {
        for (size_t sample = 0; barcodes.size() > sample; ++sample)
        {
            const std::string barcode = barcodes[sample].substr(0, barcodes[sample].length() -1);
            os << (boost::format("%s,%i,Sample%02i,artificial-genome,%s,Sample %i,N,None,Nobody,SomeProject") %
                   flowcell_ % lane % sample % barcode % sample).str() << std::endl;
        }
    }
}
*/
/*****************************************************************************/

void RunFolderGenerator::generateMatrix() const
{
    create_directory(baseCallsPath_ / "Matrix");
}

/*****************************************************************************/

void RunFolderGenerator::generatePhasing() const
{
    create_directory(baseCallsPath_ / "Phasing");
}

/*****************************************************************************/


} // namespace main
} // namespace eagle
