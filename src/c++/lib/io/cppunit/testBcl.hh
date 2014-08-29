/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#ifndef EAGLE_MODEL_TEST_BCL_HH
#define EAGLE_MODEL_TEST_BCL_HH

#include <cppunit/extensions/HelperMacros.h>

#include "io/Bcl.hh"


class TestBcl : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestBcl );
    CPPUNIT_TEST( testBclTile );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testBclTile();
};

#endif //EAGLE_MODEL_TEST_BCL_HH
