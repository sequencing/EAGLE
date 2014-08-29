################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file cxxExecutable.cmake
##
## CMake configuration file for all the c++ executables
##
## author Mauricio Varea
##
################################################################################

include (${EAGLE_GLOBALS_CMAKE})

get_filename_component(EAGLE_CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message (STATUS "Adding the c++    program subdirectory: ${EAGLE_CURRENT_DIR_NAME}")
include_directories (${EAGLE_CXX_ALL_INCLUDES})
include_directories (${CMAKE_CURRENT_BINARY_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${EAGLE_CXX_CONFIG_H_DIR})

##
## probably useless (originally defined for eland)
##
add_definitions(-D_REENTRANT)
