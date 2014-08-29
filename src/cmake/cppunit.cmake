################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file cppunit.cmake
##
## Configuration file for the cppunit subfolders
##
## author Mauricio Varea
##
################################################################################

##
## the location of the cppunit shared libraries will be needed for
## LD_LIBRARY_PATH
##

if (HAVE_CPPUNIT)
    get_filename_component(CPPUNIT_LOCATION ${HAVE_CPPUNIT} PATH)
    include_directories(${CPPUNIT_INCLUDE_DIR})
else (HAVE_CPPUNIT)
    message(FATAL_ERROR "cppunit not found")
endif (HAVE_CPPUNIT)

##
## find all the source files
##

file (GLOB EAGLE_TEST_SOURCE_LIST "*.cpp")

##
## create the targets to build the tests
##

set(EAGLE_CPPUNIT_TEST_NAME cppunitTest)
add_executable(${EAGLE_UNIQUE_PREFIX}${EAGLE_CPPUNIT_TEST_NAME} ${EAGLE_TEST_SOURCE_LIST})
set_target_properties(${EAGLE_UNIQUE_PREFIX}${EAGLE_CPPUNIT_TEST_NAME} PROPERTIES OUTPUT_NAME ${EAGLE_CPPUNIT_TEST_NAME})
add_test(cppunit_${EAGLE_UNIQUE_PREFIX}${EAGLE_CPPUNIT_TEST_NAME} "${CMAKE_CURRENT_BINARY_DIR}/${EAGLE_CPPUNIT_TEST_NAME}")

include_directories   (${EAGLE_CXX_ALL_INCLUDES} ${EAGLE_OPT_INC} "${CMAKE_SOURCE_DIR}/c++/unittest")

set(EAGLE_LINK_LIBRARIES "-lpthread")
if    (HAVE_ZLIB)
    set(EAGLE_LINK_LIBRARIES "${EAGLE_LINK_LIBRARIES} -lz")
endif (HAVE_ZLIB)
if    (HAVE_BZLIB)
    set(EAGLE_LINK_LIBRARIES "${EAGLE_LINK_LIBRARIES} -lbz2")
endif (HAVE_BZLIB)
if    (NOT EAGLE_FORCE_STATIC_LINK)
    set(EAGLE_LINK_LIBRARIES "${EAGLE_LINK_LIBRARIES} -ldl")
endif (NOT EAGLE_FORCE_STATIC_LINK)
target_link_libraries (${EAGLE_UNIQUE_PREFIX}${EAGLE_CPPUNIT_TEST_NAME}
                       eagle_${EAGLE_LIB_DIR} ${EAGLE_AVAILABLE_LIBRARIES} 
                       eagle_cppunit ${Boost_LIBRARIES} ${HAVE_CPPUNIT} 
                       ${EAGLE_LINK_LIBRARIES} ${LIBXML2_LIBRARIES} ${CPPUNIT_LIBRARY})

##
## Run some sanity check on the source file
##
foreach(EAGLE_CPPUNIT_SOURCE_FILE ${EAGLE_TEST_SOURCE_LIST})
    get_filename_component(FILE_NAME ${EAGLE_CPPUNIT_SOURCE_FILE} NAME)
    set(EAGLE_CPPUNIT_BINARY_FILE "${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME}")
    add_custom_command(TARGET ${EAGLE_UNIQUE_PREFIX}${EAGLE_CPPUNIT_TEST_NAME} 
                       PRE_BUILD
                       COMMAND ${CMAKE_SOURCE_DIR}/c++/unittest/check-source.sh ARGS ${EAGLE_CPPUNIT_BINARY_FILE}.checked ${EAGLE_CPPUNIT_SOURCE_FILE}
                       COMMENT "Sanity check on ${EAGLE_CPPUNIT_SOURCE_FILE}")
endforeach(EAGLE_CPPUNIT_SOURCE_FILE)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/RegistryNames.txt 
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/RegistryNames.txt ${CMAKE_CURRENT_BINARY_DIR}
                   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/RegistryNames.txt
                   COMMENT "Copying RegistryNames.txt for ${EAGLE_UNIQUE_PREFIX}${EAGLE_CPPUNIT_TEST_NAME}")

##
## create the targets to run the tests
##
add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/${EAGLE_CPPUNIT_TEST_NAME}.passed 
		   COMMAND export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:${CPPUNIT_LOCATION} && ./${EAGLE_CPPUNIT_TEST_NAME}
	           COMMAND touch ${EAGLE_CPPUNIT_TEST_NAME}.passed  
		   DEPENDS ${EAGLE_UNIQUE_PREFIX}${EAGLE_CPPUNIT_TEST_NAME} ${CMAKE_CURRENT_BINARY_DIR}/RegistryNames.txt
		   COMMENT "Running unit tests ${EAGLE_UNIQUE_PREFIX}${EAGLE_CPPUNIT_TEST_NAME}")
add_custom_target(${EAGLE_UNIQUE_PREFIX}passed ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${EAGLE_CPPUNIT_TEST_NAME}.passed)

##
## Copy the data directory from the source tree if available
##

find_path(${CMAKE_CURRENT_SOURCE_DIR}_DATA_DIR data PATHS ${CMAKE_CURRENT_SOURCE_DIR} NO_DEFAULT_PATH)
if (${CMAKE_CURRENT_SOURCE_DIR}_DATA_DIR)
message (STATUS "Adding the data subdirectory for the cppunits under ${EAGLE_LIB_DIR}")
    add_subdirectory (data)
    add_dependencies(${EAGLE_UNIQUE_PREFIX}passed ${EAGLE_UNIQUE_PREFIX}data)
endif (${CMAKE_CURRENT_SOURCE_DIR}_DATA_DIR)
