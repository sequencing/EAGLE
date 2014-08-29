/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#ifndef EAGLE_GENOME_TEST_VARIANT_LIST_HH
#define EAGLE_GENOME_TEST_VARIANT_LIST_HH

#include <cppunit/extensions/HelperMacros.h>

#include "genome/VariantList.hh"


class TestVariantList : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestVariantList );
    CPPUNIT_TEST( testAccess );
    CPPUNIT_TEST( testConnectivity );
    CPPUNIT_TEST_SUITE_END();
private:
    eagle::genome::VariantList VL;
public:
    TestVariantList()
    : VL( 2 )  // diploid
    {}
    void setUp();
    void tearDown();
    void testAccess();
    void testConnectivity();
};


#endif // EAGLE_GENOME_TEST_VARIANT_LIST_HH

