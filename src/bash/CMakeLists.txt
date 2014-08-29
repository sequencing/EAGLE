################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for the bash directory
##
## author Mauricio Varea
##
################################################################################

include ("${EAGLE_GLOBALS_CMAKE}")
include ("${EAGLE_MACROS_CMAKE}")
configure_files ("${CMAKE_CURRENT_SOURCE_DIR}/test" 
                 "${CMAKE_BINARY_DIR}/Testing/scripts" 
                 "*")
configure_file("test/common.sh" "${CMAKE_BINARY_DIR}/Testing/scripts/common.sh" COPYONLY)

##
## build all the internal applications for the project
##
add_subdirectory (libexec)
