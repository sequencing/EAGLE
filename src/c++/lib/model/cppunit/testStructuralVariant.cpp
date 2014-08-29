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
#include "testStructuralVariant.hh"

using eagle::model::StructuralVariant;
using eagle::model::ComplexRearrangement;
using eagle::model::Breakend;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestStructuralVariant, registryName("StructuralVariant"));


void TestStructuralVariant::setUp()
{
}

void TestStructuralVariant::tearDown()
{
}


void TestStructuralVariant::testFlags()
{
    StructuralVariant sv1("chr10", 3000UL, "A", "C");
    CPPUNIT_ASSERT(   sv1.hasSNP() );
    StructuralVariant sv2("chr10", 3000UL, "A", "ACGT");
    CPPUNIT_ASSERT(   sv2.hasInsertion() );
    StructuralVariant sv3("chr10", 3000UL, "ACGT", "A");
    CPPUNIT_ASSERT(   sv3.hasDeletion() );
}

void TestStructuralVariant::testVariant()
{   // use the unused ones for new variant types
    Breakend a1( "chr1:1000", eagle::model::Direction::NONE, "A" );
    //Breakend a2( "chr2:2000", eagle::model::Direction::NONE, "A" );
    Breakend c1( "chr1:1000", eagle::model::Direction::NONE, "C" );
    //Breakend c2( "chr2:3000", eagle::model::Direction::NONE, "C" );
    //Breakend g1( "chr1:4000", eagle::model::Direction::NONE, "G" );
    Breakend g2( "chr2:5000", eagle::model::Direction::NONE, "G" );
    Breakend t1( "chr1:6000", eagle::model::Direction::FWD, "T" );
    Breakend t2( "chr2:7000", eagle::model::Direction::REV, "T" );

    StructuralVariant sv1( ComplexRearrangement(a1,c1), eagle::model::variant::SNP );
    StructuralVariant sv2( ComplexRearrangement(g2,g2,"ACT"), eagle::model::variant::INS );
    StructuralVariant sv3( ComplexRearrangement(t1,t2), eagle::model::variant::Translocation );

    CPPUNIT_ASSERT_EQUAL( sv1.getVariant().adjacency.first.pos(), sv1.getVariant().adjacency.second.pos() );       // SNP @(pos 1000)
    CPPUNIT_ASSERT_EQUAL( string("ACT"), vector2string(sv2.getVariant().sequence) );                               // Insertion
    CPPUNIT_ASSERT_EQUAL( string(">"), sv3.getVariant().adjacency.first.dir.str() );    // Translocation (incoming == FWD)
    CPPUNIT_ASSERT_EQUAL( string("<"), sv3.getVariant().adjacency.second.dir.str() );   // Translocation (outgoing == REV)

}

