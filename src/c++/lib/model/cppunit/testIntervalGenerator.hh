/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#ifndef EAGLE_MODEL_TEST_INTERVAL_GENERATOR_HH
#define EAGLE_MODEL_TEST_INTERVAL_GENERATOR_HH

#include <cppunit/extensions/HelperMacros.h>

#include "model/IntervalGenerator.hh"


class TestIntervalGenerator : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestIntervalGenerator );
    CPPUNIT_TEST( testReproducibleRandomness );
    CPPUNIT_TEST( testLargeRandomPositions );
    CPPUNIT_TEST( testRandomIntervals );
    CPPUNIT_TEST( testUniformIntervals );
    CPPUNIT_TEST( testRandomIntervalEnd );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testReproducibleRandomness();
    void testLargeRandomPositions();
    void testRandomIntervals();
    void testUniformIntervals();
    void testRandomIntervalEnd();
};

#endif //EAGLE_MODEL_TEST_INTERVAL_GENERATOR_HH
