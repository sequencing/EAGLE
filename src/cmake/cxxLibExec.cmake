################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file cxxLibExec.cmake
##
## CMake configuration file for all the c++ libexec's
##
## author Mauricio Varea
##
################################################################################
include_directories (${EAGLE_CXX_ALL_INCLUDES})
include_directories (${CMAKE_CURRENT_BINARY_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${EAGLE_CXX_CONFIG_H_DIR})

get_filename_component(EAGLE_CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message (STATUS "Adding the c++    utility: ${EAGLE_CURRENT_DIR_NAME}")

##
## build the utility's library
##

file(GLOB EAGLE_LIBEXEC_SOURCES *.cpp *.c)
#list(REMOVE_ITEM EAGLE_LIBEXEC_SOURCES "${EAGLE_CURRENT_DIR_NAME}.cpp")

foreach (SOURCE_FILE ${EAGLE_LIBEXEC_SOURCES})
    get_filename_component(SOURCE_NAME ${SOURCE_FILE} NAME_WE)
    if (${SOURCE_NAME}_COMPILE_FLAGS)
        set_source_files_properties(${SOURCE_FILE} PROPERTIES COMPILE_FLAGS ${${SOURCE_NAME}_COMPILE_FLAGS})
    endif (${SOURCE_NAME}_COMPILE_FLAGS)
endforeach (SOURCE_FILE)

include_directories (${EAGLE_COMMON_INCLUDE})
add_executable(${EAGLE_CURRENT_DIR_NAME} ${EAGLE_LIBEXEC_SOURCES})
target_link_libraries (${EAGLE_CURRENT_DIR_NAME} ${EAGLE_AVAILABLE_LIBRARIES}
                       ${BAM_LIBRARY} ${Boost_LIBRARIES} ${LIBXML2_LIBRARIES}
                       ${EAGLE_ADDITIONAL_LIB} )
install(TARGETS ${EAGLE_CURRENT_DIR_NAME} RUNTIME DESTINATION ${EAGLE_LIBEXECDIR})


##
## build the unit tests if any (this should be mandatory really)
##

if (HAVE_CPPUNIT AND EAGLE_UNIT_TESTS)
    find_path(${CMAKE_CURRENT_SOURCE_DIR}_CPPUNIT_DIR cppunit PATHS ${CMAKE_CURRENT_SOURCE_DIR} NO_DEFAULT_PATH)
    if (${CMAKE_CURRENT_SOURCE_DIR}_CPPUNIT_DIR)
        message (STATUS "Adding the cppunit subdirectory for ${EAGLE_CURRENT_DIR_NAME}")
        add_subdirectory (cppunit)
    endif (${CMAKE_CURRENT_SOURCE_DIR}_CPPUNIT_DIR)
endif(HAVE_CPPUNIT AND EAGLE_UNIT_TESTS)

