/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#ifndef EAGLE_MODEL_TEST_PHRED_HH
#define EAGLE_MODEL_TEST_PHRED_HH

#include <cppunit/extensions/HelperMacros.h>

#include "model/Phred.hh"


class TestPhred : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestPhred );
    CPPUNIT_TEST( testQualToProb );
    CPPUNIT_TEST( testProbToQual );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testQualToProb();
    void testProbToQual();
};

#endif //EAGLE_MODEL_TEST_PHRED_HH
