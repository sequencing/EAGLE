/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#ifndef EAGLE_COMMON_TEST_EXCEPTIONS_HH
#define EAGLE_COMMON_TEST_EXCEPTIONS_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "common/Exceptions.hh"

class TestExceptions : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestExceptions );
    CPPUNIT_TEST( testErrorNumber );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testErrorNumber();
};

#endif // EAGLE_COMMON_TEST_EXCEPTIONS_HH
