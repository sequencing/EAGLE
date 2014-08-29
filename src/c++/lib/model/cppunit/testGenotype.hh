/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#ifndef EAGLE_MODEL_TEST_GENOTYPE_HH
#define EAGLE_MODEL_TEST_GENOTYPE_HH

#include <cppunit/extensions/HelperMacros.h>

#include "model/Genotype.hh"


class TestGenotype : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestGenotype );
    CPPUNIT_TEST( testPureZygosity );
    CPPUNIT_TEST( testImpureZygosity );
    CPPUNIT_TEST( testVirtualPloidy );
    CPPUNIT_TEST( testRealPloidy );
    CPPUNIT_TEST( testStreaming );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testPureZygosity();
    void testImpureZygosity();
    void testVirtualPloidy();
    void testRealPloidy();
    void testStreaming();
};

#endif // EAGLE_MODEL_TEST_GENOTYPE_HH
