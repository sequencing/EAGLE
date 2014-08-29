################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for the c++ subfolder
##
## author Mauricio Varea
##
################################################################################

set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set (EAGLE_CXX_EXECUTABLE_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxExecutable.cmake")
set (EAGLE_CXX_LIBRARY_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxLibrary.cmake")
set (EAGLE_CXX_LIBEXEC_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxLibExec.cmake")
set (EAGLE_CXX_CONFIGURE_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cxxConfigure.cmake")

set(EAGLE_CXX_CONFIG_H_DIR ${CMAKE_CURRENT_BINARY_DIR}/common)
include ("${EAGLE_CXX_CONFIGURE_CMAKE}")

if (HAVE_CPPUNIT AND EAGLE_UNIT_TESTS)
    set (EAGLE_CPPUNIT_CMAKE "${CMAKE_SOURCE_DIR}/cmake/cppunit.cmake")
    add_subdirectory (unittest)
endif (HAVE_CPPUNIT AND EAGLE_UNIT_TESTS)


##
## The include directories
##
set (EAGLE_CXX_ALL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/include")

##
## Build all the libraries for the project
##
add_subdirectory (lib)

##
## build all the applications for the project
##
add_subdirectory (bin)

##
## build all the internal applications for the project
##
add_subdirectory (libexec)

##
## build the documentation when available
##
include  (FindDoxygen)
message (STATUS "Doxygen: ${DOXYGEN_EXECUTABLE}. Dot: ${DOXYGEN_DOT_EXECUTABLE}.")
if (DOXYGEN_FOUND)
    set (DOXYFILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    message (STATUS "Creating Doxygen config file: ${DOXYFILE}")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${DOXYFILE} @ONLY IMMEDIATE)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOXYFILE})
endif (DOXYGEN_FOUND)
