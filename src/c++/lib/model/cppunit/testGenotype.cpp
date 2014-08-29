/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#include <iostream>
#include <string>
#include <sstream>
#include <boost/assign.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testGenotype.hh"

using eagle::model::Genotype;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestGenotype, registryName("Genotype"));

void TestGenotype::setUp()
{
}

void TestGenotype::tearDown()
{
}


void TestGenotype::testPureZygosity()
{
    Genotype gt;
    CPPUNIT_ASSERT(  gt.isHomozygousRef() );
    gt.set(0);
    CPPUNIT_ASSERT(  gt.isHomozygousRef() );
    CPPUNIT_ASSERT(! gt.isHomozygousDiff() );
    CPPUNIT_ASSERT(  gt.isHomozygous() );
    CPPUNIT_ASSERT(! gt.isHeterozygous() );
    gt.setPloidy(2);
    gt.set(1);
    CPPUNIT_ASSERT(! gt.isHomozygousRef() );
    CPPUNIT_ASSERT(! gt.isHomozygousDiff() );
    CPPUNIT_ASSERT(! gt.isHomozygous() );
    CPPUNIT_ASSERT(! gt.isHeterozygous() );
    gt.set(3);
    gt.reset(0);
    CPPUNIT_ASSERT(! gt.isHomozygousRef() );
    CPPUNIT_ASSERT(  gt.isHomozygousDiff() );
    CPPUNIT_ASSERT(  gt.isHomozygous() );
    CPPUNIT_ASSERT(! gt.isHeterozygous() );
    gt.reset(1);
    CPPUNIT_ASSERT(! gt.isHomozygousRef() );
    CPPUNIT_ASSERT(! gt.isHomozygousDiff() );
    CPPUNIT_ASSERT(! gt.isHomozygous() );
    CPPUNIT_ASSERT(  gt.isHeterozygous() );
}

void TestGenotype::testImpureZygosity()
{
    Genotype gt(2);
    gt.set(0);
    gt.set(1);
    gt.set(-1);
    CPPUNIT_ASSERT(! gt.isHomozygousRef() );
    CPPUNIT_ASSERT(! gt.isHomozygousDiff() );
    CPPUNIT_ASSERT(! gt.isHomozygous() );
    CPPUNIT_ASSERT(! gt.isHeterozygous() );
}

void TestGenotype::testVirtualPloidy()
{
    Genotype gt;
    CPPUNIT_ASSERT_EQUAL( 1u, gt.minPloidy() );
    CPPUNIT_ASSERT_EQUAL( 1u, gt.maxPloidy() );
    gt.setPloidy(4);
    gt.set(0);
    gt.set(-1);
    CPPUNIT_ASSERT_EQUAL( 1u, gt.minPloidy() );
    CPPUNIT_ASSERT_EQUAL( 4u, gt.maxPloidy() );
    gt.set(5);
    CPPUNIT_ASSERT_EQUAL( 5u, gt.minPloidy() );
    CPPUNIT_ASSERT_EQUAL( 5u, gt.maxPloidy() );
}

void TestGenotype::testRealPloidy()
{
    Genotype gt(3);
    gt.set(2);
    gt.set(-1);
    gt.set(7);
    CPPUNIT_ASSERT_EQUAL( 2u, gt.minPloidy() );
    CPPUNIT_ASSERT_EQUAL( 3u, gt.getPloidy() );
    CPPUNIT_ASSERT_EQUAL( 7u, gt.maxPloidy() );
    CPPUNIT_ASSERT_EQUAL( 2u, gt.altSize() );
}

void TestGenotype::testStreaming()
{
    stringstream ss1;
    Genotype gt(5);
    gt.set(-1);
    gt.set(0);
    ss1 << gt;
    CPPUNIT_ASSERT_EQUAL( string("0/0/0/0/0"), ss1.str() );

    stringstream ss2;
    gt.set(2);
    gt.set(7);
    gt.set(3);
    ss2 << gt;
    //CPPUNIT_ASSERT_EQUAL( string("1/1/0/0/0/1"), ss2.str() );
    CPPUNIT_ASSERT_EQUAL( string("0/1/1/0/0/0/1"), ss2.str() );
}
