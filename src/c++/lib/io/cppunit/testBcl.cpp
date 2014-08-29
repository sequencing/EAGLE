/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#include <iostream>
#include <string>
#include <vector>
#include <boost/assign.hpp>

using namespace std;
using boost::assign::list_of;

#include "Helpers.hh"

#include "RegistryName.hh"
#include "testBcl.hh"

using eagle::io::BclTile;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestBcl, registryName("Bcl"));


void TestBcl::setUp()
{
}

void TestBcl::tearDown()
{
}


void TestBcl::testBclTile()
{
    // Create a tile
    const unsigned int expectedReadCount = 5;
    const unsigned int clusterLength     = 101;
    const std::string filenameTemplate      = "";
    const std::string statsFilenameTemplate = "";
    const std::string filterFilename        = "";
    const std::string posFilename           = "";
    const std::string controlFilename       = "";

    BclTile tile( expectedReadCount, clusterLength, filenameTemplate, statsFilenameTemplate, filterFilename, posFilename, controlFilename, false );

    // Test that we can fill it
    const char bufCluster[102] = "12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901";
    CPPUNIT_ASSERT_NO_THROW( tile.addClusterToRandomLocation( bufCluster ) );
    CPPUNIT_ASSERT_NO_THROW( tile.addClusterToRandomLocation( bufCluster ) );
    CPPUNIT_ASSERT_NO_THROW( tile.addClusterToRandomLocation( bufCluster ) );
    CPPUNIT_ASSERT_NO_THROW( tile.addClusterToRandomLocation( bufCluster ) );
    CPPUNIT_ASSERT_NO_THROW( tile.addClusterToRandomLocation( bufCluster ) ); // expectedReadCount == 5 => tile is now full

    // Test that we cannot overfill it
    CPPUNIT_ASSERT_THROW( tile.addClusterToRandomLocation( bufCluster ), eagle::common::OutOfLimitsException ); // too many clusters added to tile
}
