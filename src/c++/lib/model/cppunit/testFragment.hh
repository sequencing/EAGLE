/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#ifndef EAGLE_MODEL_TEST_FRAGMENT_HH
#define EAGLE_MODEL_TEST_FRAGMENT_HH

#include <cppunit/extensions/HelperMacros.h>

#include "model/Fragment.hh"


class TestFragment : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestFragment );
    CPPUNIT_TEST( testFragment );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testFragment();
};

#endif //EAGLE_MODEL_TEST_FRAGMENT_HH
