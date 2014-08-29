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
#include "testIntervalGenerator.hh"

using eagle::model::RandomSequenceGenerator;
using eagle::model::RandomIntervalGenerator;
using eagle::model::UniformIntervalGenerator;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestIntervalGenerator, registryName("IntervalGenerator"));


void TestIntervalGenerator::setUp()
{
}

void TestIntervalGenerator::tearDown()
{
}


void TestIntervalGenerator::testReproducibleRandomness()
{
    // This test checks if we can expect reproducible results from the random number generator
    RandomSequenceGenerator gen( 5, 100 );

    CPPUNIT_ASSERT_EQUAL( gen.getNext(), 3lu );
    CPPUNIT_ASSERT_EQUAL( gen.getNext(), 23lu );
    CPPUNIT_ASSERT_EQUAL( gen.getNext(), 29lu );
    CPPUNIT_ASSERT_EQUAL( gen.getNext(), 36lu );
    CPPUNIT_ASSERT_EQUAL( gen.getNext(), 42lu );
    CPPUNIT_ASSERT_THROW( gen.getNext(), eagle::common::EagleException );
}

void TestIntervalGenerator::testLargeRandomPositions()
{
    // This test relies on reproducible results from the random number generator (previous test)
    // to test that large numbers (larger than 32bits) don't get truncated anywhere
    RandomSequenceGenerator gen( 3, 15000000000000000000ul ); // large unsigned long, 64 bits or more

    CPPUNIT_ASSERT_EQUAL( gen.getNext(),  6263893078224449536ul );
    CPPUNIT_ASSERT_EQUAL( gen.getNext(),  9941931729156847616ul );
    CPPUNIT_ASSERT_EQUAL( gen.getNext(), 11114242261213317120ul );
}

void TestIntervalGenerator::testRandomIntervals()
{
    // Check that it works with small chromosomes in the middle of large ones
    const std::vector<unsigned long> contigLengths        = list_of(10000000000ul)(50ul)(100ul)(10000000000ul)(50ul);
    const unsigned long              readCount            = 1000;
    const unsigned int               minFragmentLength    = 10;
    const unsigned int               medianFragmentLength = 20;
    const unsigned int               maxFragmentLength    = 30;
    RandomIntervalGenerator gen( contigLengths, readCount, minFragmentLength, medianFragmentLength, maxFragmentLength, false);

    // Check that numbers get generated in increasing order
    unsigned long lastPos = 0;
    for (unsigned int i=0; i<readCount; ++i)
    {
        unsigned long pos = gen.getNext().first;
        CPPUNIT_ASSERT( pos >= lastPos );
        lastPos = pos;
    }
}

void TestIntervalGenerator::testUniformIntervals()
{
    // Check that it works with small chromosomes in the middle of large ones
    const std::vector<unsigned long> contigLengths        = list_of(14ul)(100000000000ul)(50ul)(100ul)(10000000000ul)(50ul);
    const unsigned int               medianFragmentLength = 10;
    const double                     step                 = 0.75;
    unsigned long                    readCount            = 1000;
    UniformIntervalGenerator gen( contigLengths, medianFragmentLength, step, readCount, false);

    CPPUNIT_ASSERT_EQUAL( readCount, 146666666880ul );
    CPPUNIT_ASSERT_EQUAL( gen.getNext().first,  0ul );
    CPPUNIT_ASSERT_EQUAL( gen.getNext().first,  0ul ); // Here and below, we check that it takes into account the floating point step
    CPPUNIT_ASSERT_EQUAL( gen.getNext().first,  1ul );
    CPPUNIT_ASSERT_EQUAL( gen.getNext().first,  2ul );
    CPPUNIT_ASSERT_EQUAL( gen.getNext().first,  3ul );
    CPPUNIT_ASSERT_EQUAL( gen.getNext().first,  3ul );
    CPPUNIT_ASSERT_EQUAL( gen.getNext().first,  4ul );
    CPPUNIT_ASSERT_EQUAL( gen.getNext().first, 14ul ); // Here we check that it jumps to next contig because of the fragment length
    CPPUNIT_ASSERT_EQUAL( gen.getNext().first, 14ul ); // Here we check that the jump did reset the position to 14.0 and not to 14.25
    CPPUNIT_ASSERT_EQUAL( gen.getNext().first, 15ul );
    CPPUNIT_ASSERT_EQUAL( gen.getNext().first, 16ul );
    CPPUNIT_ASSERT_EQUAL( gen.getNext().first, 17ul );
}

void TestIntervalGenerator::testRandomIntervalEnd()
{
    // We want to check that random intervals get properly enumerated and that the conversion enumeration to interval is not bugged
    // The difficulty is that these intervals are kept private, 
    // and that the conversion only happens when an interval gets generated, which is done in a random way
    // This tests tries to go through all the random intervals by building an environment with only a few intervals
    // We generate many intervals to achieve a complete coverage, and we check that they are all present in the expected order
    const std::vector<unsigned long> contigLengths        = list_of(6ul)(4ul);
    const unsigned long              readCount            = 1000;
    const unsigned int               minFragmentLength    = 2;
    const unsigned int               medianFragmentLength = 3;
    const unsigned int               maxFragmentLength    = 4;

    const std::vector< std::pair< unsigned long, unsigned int > > expectedIntervals = list_of
        (std::make_pair(0ul,2))
        (std::make_pair(0ul,3))
        (std::make_pair(0ul,4))

        (std::make_pair(1ul,2))
        (std::make_pair(1ul,3))
        (std::make_pair(1ul,4))

        (std::make_pair(2ul,2))
        (std::make_pair(2ul,3))
        (std::make_pair(2ul,4))

        (std::make_pair(3ul,2))
        (std::make_pair(3ul,3))

        (std::make_pair(4ul,2))


        (std::make_pair(6ul,2))
        (std::make_pair(6ul,3))
        (std::make_pair(6ul,4))

        (std::make_pair(7ul,2))
        (std::make_pair(7ul,3))

        (std::make_pair(8ul,2))
         ;

    RandomIntervalGenerator gen( contigLengths, readCount, minFragmentLength, medianFragmentLength, maxFragmentLength, false);

    // Check that numbers get generated in increasing order
    unsigned int lastIntervalNum = 0;
    for (unsigned int i=0; i<readCount; ++i)
    {
        std::pair< unsigned long, unsigned int > interval = gen.getNext();
        //        std::clog << "int: " << interval.first << "," << interval.second << std::endl;
        if (interval != expectedIntervals[lastIntervalNum])
        {
            if (lastIntervalNum < expectedIntervals.size()-1)
            {
                ++lastIntervalNum;
            }
            if (interval != expectedIntervals[lastIntervalNum])
            {
                std::cerr << "expected  interval: [" << expectedIntervals[lastIntervalNum].first << "," << expectedIntervals[lastIntervalNum].second << "]" << std::endl;
                std::cerr << "generated interval: [" << interval.first << "," << interval.second << "]" << std::endl;
                CPPUNIT_ASSERT( false );
            }
        }
    }
}

