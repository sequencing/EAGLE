################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for the libexec subfolder
##
## author Mauricio Varea
##
################################################################################

##
## List of EAGLE utilities
##
file (GLOB EAGLE_ALL_UTILITIES [a-z][a-zA-Z0-9]*)

##
## Build all the utilities for the project
## EAGLE_AVAILABLE_UTILITIES is incrementally updated
##

set (EAGLE_AVAILABLE_UTILITIES "")
set (EAGLE_ALL_UTILITIES_DIRS "")
foreach (EAGLE_LIBEXEC ${EAGLE_ALL_UTILITIES})
    add_subdirectory(${EAGLE_LIBEXEC})
    set(EAGLE_AVAILABLE_UTILITIES ${EAGLE_LIBEXEC} ${EAGLE_AVAILABLE_UTILITIES} )
endforeach (EAGLE_LIBEXEC)

set (EAGLE_AVAILABLE_UTILITIES ${EAGLE_AVAILABLE_UTILITIES} PARENT_SCOPE)
