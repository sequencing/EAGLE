################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for the share subdirectory
##
## author Mauricio Varea
##
################################################################################

#add_subdirectory (verifyBoost)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../CTestCustom.cmake" "${CMAKE_BINARY_DIR}/" COPYONLY)