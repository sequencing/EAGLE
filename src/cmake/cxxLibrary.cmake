################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file cxxLibrary.cmake
##
## CMake configuration file for all the c++ libraries
##
## author Mauricio Varea
##
################################################################################
include_directories (${EAGLE_CXX_ALL_INCLUDES})
include_directories (${CMAKE_CURRENT_BINARY_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${EAGLE_CXX_CONFIG_H_DIR})

get_filename_component(EAGLE_CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message (STATUS "Adding the c++    library subdirectory: ${EAGLE_CURRENT_DIR_NAME}")

##
## Some generators (VS) require all targets to be unique across the project.
## Therefore, a unique prefix is needed to create the target names which are
## shared across libraries
##

string(REGEX REPLACE ${CMAKE_SOURCE_DIR}/c[+][+]/ "" TMP1 ${CMAKE_CURRENT_SOURCE_DIR}/)
string(REGEX REPLACE "/" "_" EAGLE_UNIQUE_PREFIX ${TMP1})

##
## build the library
##

file(GLOB EAGLE_LIBRARY_SOURCES *.cpp *.c)
foreach (SOURCE_FILE ${EAGLE_LIBRARY_SOURCES})
    get_filename_component(SOURCE_NAME ${SOURCE_FILE} NAME_WE)
    if (${SOURCE_NAME}_COMPILE_FLAGS)
        set_source_files_properties(${SOURCE_FILE} PROPERTIES COMPILE_FLAGS ${${SOURCE_NAME}_COMPILE_FLAGS})
    endif (${SOURCE_NAME}_COMPILE_FLAGS)
endforeach (SOURCE_FILE)

include_directories (${EAGLE_COMMON_INCLUDE} )
add_library         (eagle_${EAGLE_LIB_DIR} STATIC ${EAGLE_LIBRARY_SOURCES})

##
## build the unit tests if any (this should be mandatory really)
##

if (HAVE_CPPUNIT AND EAGLE_UNIT_TESTS)
    find_path(${CMAKE_CURRENT_SOURCE_DIR}_CPPUNIT_DIR cppunit PATHS ${CMAKE_CURRENT_SOURCE_DIR} NO_DEFAULT_PATH)
    if (${CMAKE_CURRENT_SOURCE_DIR}_CPPUNIT_DIR)
        message (STATUS "Adding the cppunit subdirectory for ${EAGLE_LIB_DIR}")
        add_subdirectory (cppunit)
    endif (${CMAKE_CURRENT_SOURCE_DIR}_CPPUNIT_DIR)
endif(HAVE_CPPUNIT AND EAGLE_UNIT_TESTS)

