/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#ifndef EAGLE_MODEL_TEST_VCF_HH
#define EAGLE_MODEL_TEST_VCF_HH

#include <cppunit/extensions/HelperMacros.h>

#include "io/Vcf.hh"


class TestVcf : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestVcf );
    CPPUNIT_TEST( testVcfParserSnp );
    CPPUNIT_TEST( testVcfParserInsertions );
    CPPUNIT_TEST( testVcfParserDeletions );
    CPPUNIT_TEST( testVcfParserIndels );
    CPPUNIT_TEST( testVcfParserTranslocations );
    CPPUNIT_TEST( testVcfParserInversions );
    CPPUNIT_TEST( testVcfParserDuplications );
    CPPUNIT_TEST( testVcfParserQualityField );
    CPPUNIT_TEST( testVcfParserFilterField );
    CPPUNIT_TEST( testVcfParserInfoField );
    CPPUNIT_TEST_SUITE_END();
public:
    void setUp();
    void tearDown();
    void testVcfParserSnp();
    void testVcfParserInsertions();
    void testVcfParserDeletions();
    void testVcfParserIndels();
    void testVcfParserTranslocations();
    void testVcfParserInversions();
    void testVcfParserDuplications();
    void testVcfParserQualityField();
    void testVcfParserFilterField();
    void testVcfParserInfoField();
};

#endif //EAGLE_MODEL_TEST_VCF_HH
