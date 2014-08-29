/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#include <iostream>
#include <string>
#include <vector>
#include <boost/assign.hpp>

using namespace std;
using boost::assign::list_of;

#include "Helpers.hh"

#include "RegistryName.hh"
#include "testPhred.hh"

using eagle::model::Phred;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestPhred, registryName("Phred"));


void TestPhred::setUp()
{
}

void TestPhred::tearDown()
{
}

void TestPhred::testQualToProb()
{
    CPPUNIT_ASSERT_EQUAL( Phred::qualToProb(0), 1.0 );
    CPPUNIT_ASSERT_EQUAL( Phred::qualToProb(1), pow(10, -0.1) );
    CPPUNIT_ASSERT_EQUAL( Phred::qualToProb(10), 0.1 );
    CPPUNIT_ASSERT_EQUAL( Phred::qualToProb(20), 0.01 );
    CPPUNIT_ASSERT_EQUAL( Phred::qualToProb(30), 0.001 );
    CPPUNIT_ASSERT_EQUAL( Phred::qualToProb(40), 0.0001 );
    CPPUNIT_ASSERT_EQUAL( Phred::qualToProb(50), 0.00001 );
    CPPUNIT_ASSERT_THROW( Phred::qualToProb(51), eagle::common::EagleException );
}

void TestPhred::testProbToQual()
{
    CPPUNIT_ASSERT_EQUAL( Phred::probToQual(1.0),            0u );
    CPPUNIT_ASSERT_EQUAL( Phred::probToQual(pow(10, -0.1)),  1u );
    CPPUNIT_ASSERT_EQUAL( Phred::probToQual(0.1),           10u );
    CPPUNIT_ASSERT_EQUAL( Phred::probToQual(0.01),          20u );
    CPPUNIT_ASSERT_EQUAL( Phred::probToQual(0.001),         30u );
    CPPUNIT_ASSERT_EQUAL( Phred::probToQual(0.0001),        40u );
    CPPUNIT_ASSERT_EQUAL( Phred::probToQual(0.00001),       50u );
    CPPUNIT_ASSERT_EQUAL( Phred::probToQual(0.0),           51u );

    double between10and11 = (Phred::qualToProb(10) + Phred::qualToProb(11)) / 2;
    CPPUNIT_ASSERT_EQUAL( Phred::probToQual(between10and11), 11u );
}
