################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for the c++/unittest subfolder
##
## author Come Raczy
##
################################################################################

include_directories(${CPPUNIT_INCLUDE_DIR})
add_library(eagle_cppunit cppunitTest.cpp RegistryName.cpp)
