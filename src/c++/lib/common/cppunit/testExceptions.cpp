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
#include "testExceptions.hh"

using namespace eagle::common;


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestExceptions, registryName("Exceptions"));

void TestExceptions::setUp()
{
}

void TestExceptions::tearDown()
{
}


void TestExceptions::testErrorNumber()
{
    ExceptionData e(1337, "exception unit test");
    CPPUNIT_ASSERT_EQUAL(1337, e.getErrorNumber());
    CPPUNIT_ASSERT_EQUAL(string("exception unit test"), e.getMessage());
}

