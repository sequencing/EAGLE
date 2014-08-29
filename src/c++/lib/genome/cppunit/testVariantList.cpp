/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#include <iostream>
#include <string>
#include <boost/assign.hpp>

using namespace std;
using boost::assign::list_of;

#include "Helpers.hh"

#include "RegistryName.hh"
#include "testVariantList.hh"

using eagle::genome::VariantList;
using eagle::genome::Event;
using eagle::model::StructuralVariant;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestVariantList, registryName("VariantList"));


void TestVariantList::setUp()
{
    VL.push_back( Event( StructuralVariant("chr10", 4000UL, "A", "C") ) );
    VL.push_back( Event( StructuralVariant("chr10", 5000UL, "A", "ACGT") ) );
    VL.push_back( Event( StructuralVariant("chr10", 6000UL, "ACGT", "A") ) );
    /*
    VL.push_front( Event( StructuralVariant("chr10", 3000UL, "ACGT", "A") ) );
    VL.push_front( Event( StructuralVariant("chr10", 2000UL, "A", "ACGT") ) );
    VL.push_front( Event( StructuralVariant("chr10", 1000UL, "A", "C") ) );
    */
}

void TestVariantList::tearDown()
{
    while (!VL.empty())
    {
        VL.pop_back();
    }
}


void TestVariantList::testAccess()
{
    Event third( StructuralVariant("chr10", 3000UL, "ACGT", "A") );
    /*
    CPPUNIT_ASSERT_EQUAL(third, *(VL.before() + 2) );
    CPPUNIT_ASSERT_EQUAL(third, *(VL.last() - 3) );
    CPPUNIT_ASSERT_EQUAL(third, *(VL.afterLast() - 4) );
    */
}

void TestVariantList::testConnectivity()
{
    // 999 bases are transferred when applying deletion@3000
    //CPPUNIT_ASSERT_EQUAL(999UL, (unsigned long)(VL.to(VL.before()+2) - VL.from(VL.before()+1)) );
    // 996 bases are transferred when applying SNP@4000
    //CPPUNIT_ASSERT_EQUAL(996UL, (unsigned long)(VL.to(VL.before()+3) - VL.from(VL.before()+2)) );
}

