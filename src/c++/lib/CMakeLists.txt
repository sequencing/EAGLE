################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for the lib subfolder
##
## author Mauricio Varea
##
################################################################################

##
## List of EAGLE libraries
##
set (EAGLE_ALL_LIBRARIES common model io genome main)

##
## Build all the libraries for the project
## EAGLE_AVAILABLE_LIBRARIES is incrementally updated
##

set (EAGLE_AVAILABLE_LIBRARIES "")
set (EAGLE_ALL_LIBRARY_DIRS "")
foreach (EAGLE_LIB_DIR ${EAGLE_ALL_LIBRARIES})
    add_subdirectory(${EAGLE_LIB_DIR})
    set(EAGLE_AVAILABLE_LIBRARIES eagle_${EAGLE_LIB_DIR} ${EAGLE_AVAILABLE_LIBRARIES} )
endforeach (EAGLE_LIB_DIR)

set (EAGLE_AVAILABLE_LIBRARIES ${EAGLE_AVAILABLE_LIBRARIES} PARENT_SCOPE)

