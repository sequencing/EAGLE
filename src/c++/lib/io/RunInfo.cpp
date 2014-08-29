/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Component to read/write RunInfo.xml files.
 **
 ** \author Lilian Janin
 **/

#include <boost/format.hpp>
#include "io/RunInfo.hh"
#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

namespace eagle
{
namespace io
{

void RunInfo::parse( const boost::filesystem::path &filename )
{
    boost::property_tree::ptree runinfo_xml;
    read_xml( filename.string(), runinfo_xml );
    try {
        boost::property_tree::ptree &entryRun = runinfo_xml.get_child("RunInfo.Run");
        runId          = entryRun.get<string>("<xmlattr>.Id"            , "DEFAULT_RUN_ID");
        runNumber      = entryRun.get<string>("<xmlattr>.Number"        , "0");
        tileNameMethod = entryRun.get<string>("<xmlattr>.TileNameMethod", "1");

        flowcell       = entryRun.get<string>("Flowcell"                , "FC1234XXX");

        boost::property_tree::ptree &entryFlowcellLayout = entryRun.get_child("FlowcellLayout");
        laneCount    = entryFlowcellLayout.get<unsigned int>("<xmlattr>.LaneCount");
        surfaceCount = entryFlowcellLayout.get<unsigned int>("<xmlattr>.SurfaceCount", 1);
        swathCount   = entryFlowcellLayout.get<unsigned int>("<xmlattr>.SwathCount"  , 1);
        tileCount    = entryFlowcellLayout.get<unsigned int>("<xmlattr>.TileCount");

        boost::property_tree::ptree &entryReads = entryRun.get_child("Reads");
        BOOST_FOREACH (boost::property_tree::ptree::value_type &read, entryReads)
        {
            if (read.first == "Read")
            {
                ReadDescription readDesc;
                readDesc.firstCycle = read.second.get<unsigned int>("<xmlattr>.FirstCycle");
                readDesc.lastCycle = read.second.get<unsigned int>("<xmlattr>.LastCycle");
                boost::optional<boost::property_tree::ptree &> optionalEntryIndex = read.second.get_child_optional("Index");
                readDesc.isIndex = optionalEntryIndex.is_initialized();

                reads.push_back(readDesc);
            }
        }
    }
    catch (boost::property_tree::ptree_bad_path &e)
    {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("XML.RunInfo",
                             (boost::format("*** Could not parse %s ***") % filename).str() ));
    }

    EAGLE_DEBUG(0, "Parsed RunInfo.xml data:" );
    EAGLE_DEBUG(7, "runId=" << runId );
    EAGLE_DEBUG(7, "runNumber=" << runNumber );
    EAGLE_DEBUG(7, "tileNameMethod=" << tileNameMethod );
    EAGLE_DEBUG(7, "flowcell=" << flowcell );
    EAGLE_DEBUG(7, "laneCount=" << laneCount );
    EAGLE_DEBUG(7, "surfaceCount=" << surfaceCount );
    EAGLE_DEBUG(7, "swathCount=" << swathCount );
    EAGLE_DEBUG(7, "tileCount=" << tileCount );
    EAGLE_DEBUG(7, "reads count: " << reads.size() );
#ifdef EAGLE_DEBUG_MODE
    BOOST_FOREACH(ReadDescription &rd, reads)
    {
        EAGLE_DEBUG(14, (boost::format("read: {%d,%d,%d}") % rd.firstCycle % rd.lastCycle % rd.isIndex).str() );
    }
#endif
}

} // namespace io
} // namespace eagle
