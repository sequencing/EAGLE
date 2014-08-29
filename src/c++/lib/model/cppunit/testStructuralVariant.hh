/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#ifndef EAGLE_MODEL_TEST_STRUCTURAL_VARIANT_HH
#define EAGLE_MODEL_TEST_STRUCTURAL_VARIANT_HH

#include <cppunit/extensions/HelperMacros.h>

#include "model/StructuralVariant.hh"


class TestStructuralVariant : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestStructuralVariant );
    CPPUNIT_TEST( testFlags );
    CPPUNIT_TEST( testVariant );
    CPPUNIT_TEST_SUITE_END();
private:
public:
    void setUp();
    void tearDown();
    void testFlags();
    void testVariant();
};

#endif // EAGLE_MODEL_TEST_STRUCTURAL_VARIANT_HH

