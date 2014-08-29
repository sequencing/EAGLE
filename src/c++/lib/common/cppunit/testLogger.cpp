/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#include <iostream>
#include <sstream>
#include <string>

using namespace std;

#include "RegistryName.hh"
#include "testLogger.hh"

using namespace eagle::common;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestLogger, registryName("Logger"));

static const size_t MiliSec(1000);
static const size_t Sec(1000 * MiliSec);
static const size_t Min(60 * Sec);
static const size_t Hour(60 * Min);

void TestLogger::setUp()
{
}

void TestLogger::tearDown()
{
}


void TestLogger::testTime()
{
    CPPUNIT_ASSERT_EQUAL( string("999ms"),                  displayTime(static_cast<size_t>(999 * MiliSec)) );
    CPPUNIT_ASSERT_EQUAL( string("1000ms (0h:0m:1s)"),      displayTime(static_cast<size_t>(  1 * Sec)) );
    CPPUNIT_ASSERT_EQUAL( string("60000ms (0h:1m:0s)"),     displayTime(static_cast<size_t>(  1 * Min)) );
    CPPUNIT_ASSERT_EQUAL( string("3600000ms (1h:0m:0s)"),   displayTime(static_cast<size_t>(  1 * Hour)) );
    CPPUNIT_ASSERT_EQUAL( string("9045500ms (2h:30m:45s)"), displayTime(static_cast<size_t>(2 * Hour + 30 * Min + 45 * Sec + 500 * MiliSec)) );
}

