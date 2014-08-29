/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#ifndef EAGLE_COMMON_TEST_LOGGER_HH
#define EAGLE_COMMON_TEST_LOGGER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "common/Logger.hh"

class TestLogger : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestLogger );
    CPPUNIT_TEST( testTime );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testTime();
};

#endif // EAGLE_COMMON_TEST_LOGGER_HH
