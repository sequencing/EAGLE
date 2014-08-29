/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#ifndef EAGLE_MODEL_TEST_CONTIG_HH
#define EAGLE_MODEL_TEST_CONTIG_HH

#include <cppunit/extensions/HelperMacros.h>

#include "model/Contig.hh"


class TestContig : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestContig );
    CPPUNIT_TEST( testHeadReset );
    CPPUNIT_TEST( testHeadGet );
    CPPUNIT_TEST( testHeadCreate );
    CPPUNIT_TEST( testBodyBlock );
    //CPPUNIT_TEST( testBodySingle );
    //CPPUNIT_TEST( testBodyManipulation );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testHeadReset();
    void testHeadGet();
    void testHeadCreate();
    void testBodyBlock();
    //void testBodySingle();
    //void testBodyManipulation();
};

#endif // EAGLE_MODEL_TEST_CONTIG_HH

