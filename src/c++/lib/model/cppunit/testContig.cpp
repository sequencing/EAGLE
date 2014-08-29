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
#include "testContig.hh"

using eagle::model::Contig;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestContig, registryName("Contig"));

static const Contig e_coli("gi|49175990|ref|NC_000913.2| Escherichia coli K12, complete genome");
static const Contig human("chromosomeN");

static const vector<char> E_coli_10 = string2vector("GACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTT");

void TestContig::setUp()
{
}

void TestContig::tearDown()
{
}


void TestContig::testHeadReset()
{
    Contig dummy("my_head");
    CPPUNIT_ASSERT_EQUAL(string("my_head"), dummy.name());
    dummy.reset();
    CPPUNIT_ASSERT_EQUAL(string(""), dummy.name());
}

void TestContig::testHeadGet()
{
    CPPUNIT_ASSERT_EQUAL(string("gi|49175990|ref|NC_000913.2|"), e_coli.id());
    CPPUNIT_ASSERT_EQUAL(string("Escherichia coli K12, complete genome"), e_coli.remainder());
    CPPUNIT_ASSERT_EQUAL(string("chromosomeN"), human.id());
    CPPUNIT_ASSERT_EQUAL(string(""), human.remainder());
}

void TestContig::testHeadCreate()
{
    Contig clone1;
    Contig clone2;
    clone1.name("gi|49175990|ref|NC_000913.2|","Escherichia coli K12, complete genome");
    CPPUNIT_ASSERT_EQUAL(e_coli.name(), clone1.name());
    clone2.name("chromosomeN");
    CPPUNIT_ASSERT_EQUAL(human.name(), clone2.name());
}

void TestContig::testBodyBlock()
{
    Contig line10;
    CPPUNIT_ASSERT_EQUAL( 20UL, line10.append( subv( E_coli_10, 0, 20 )) );
    CPPUNIT_ASSERT_EQUAL( 20UL, line10.append( subv( E_coli_10, 20, 20 )) );
    CPPUNIT_ASSERT_EQUAL( string("GACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCG"),           vector2string(line10.read()) );
    CPPUNIT_ASSERT_EQUAL( 10UL, line10.append( subv( E_coli_10, 40, 10 ), true) );      // comp( CAATTGAAAA )  -- not reversed!
    //                            12345678901234567890123456789012345678901234567890
    CPPUNIT_ASSERT_EQUAL( string("GACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGGTTAACTTTT"), vector2string(line10.read()) );
    CPPUNIT_ASSERT_EQUAL(                    string("CAGCCGGGGTTCCCGCTGGCGGTTAACTTTT"), vector2string(line10.read(20)) );
    CPPUNIT_ASSERT_EQUAL(          string("CGCCGCCGCCCAGCCGGGGT"),                      vector2string(line10.read(10,30)) );
}

/*
void TestContig::testBodySingle()
{
//    test put() and get()
}

void TestContig::testBodyManipulation()
{
//    test ins() and del()
}
*/
