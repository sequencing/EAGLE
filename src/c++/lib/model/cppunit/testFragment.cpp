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
#include "testFragment.hh"

using eagle::model::Fragment;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestFragment, registryName("Fragment"));


void TestFragment::setUp()
{
}

void TestFragment::tearDown()
{
}


void TestFragment::testFragment()
{
}
