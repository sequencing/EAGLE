/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Main program used for all the cppunit tests
 **
 ** \author Come Raczy
 **/

#include <sstream>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include "RegistryName.hh"

int main()
{

  CppUnit::TextUi::TestRunner runner;
  // First add the tests from the named registries in the right order
  // To add/remove/modify a registry name, or to change its sequence
  // number, edit RegistryName.cpp
  for (std::vector<std::string>::const_iterator name = getRegistryNameList().begin();
       getRegistryNameList().end() != name; ++name)
  {
    CppUnit::Test *namedSuite = CppUnit::TestFactoryRegistry::getRegistry(*name).makeTest();
    runner.addTest(namedSuite);
  }

  // Add the top level (unnamed) suite from the list of tests to run
  CppUnit::Test *suite = CppUnit::TestFactoryRegistry::getRegistry().makeTest();
  runner.addTest( suite );

  // Change the default outputter to a compiler error format outputter
  runner.setOutputter( new CppUnit::CompilerOutputter( &runner.result(),
                                                       std::cerr ) );
  // Run the tests.
  bool wasSuccessful = runner.run();

  // Return error code 1 if the one of test failed.
  return wasSuccessful ? 0 : 1;
}
