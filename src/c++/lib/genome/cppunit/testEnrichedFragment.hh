/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#ifndef EAGLE_GENOME_TEST_ENRICHED_FRAGMENT_HH
#define EAGLE_GENOME_TEST_ENRICHED_FRAGMENT_HH

#include <cppunit/extensions/HelperMacros.h>

#include "genome/EnrichedFragment.hh"


class TestEnrichedFragment : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestEnrichedFragment );
    CPPUNIT_TEST( testEnrichedFragment );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testEnrichedFragment();
};

#endif //EAGLE_GENOME_TEST_FRAGMENT_HH
