/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Management of the registry names for the cppunit tests.
 **
 ** \author Come Raczy
 **/

#ifndef EAGLE_UNIT_TEST_REGISTRY_NAME
#define EAGLE_UNIT_TEST_REGISTRY_NAME

#include <stdexcept>
#include <string>
#include <vector>

const std::vector<std::string> &getRegistryNameList();
std::string registryName(const std::string &name) throw (std::invalid_argument);

#endif // EAGLE_UNIT_TEST_REGISTRY_NAME
