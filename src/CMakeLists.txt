################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Top level cmake configuration file
##
## author Mauricio Varea
##
################################################################################

project(EAGLE)

cmake_minimum_required(VERSION 2.8.0)
INCLUDE(CTest)

if (DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
else()
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
endif()

set(EAGLE_NAME "Enhanced Artificial Genome Engine (EAGLE)")
set(EAGLE_COPYRIGHT "Copyright (c) 2014 Illumina, Inc.")

# The EAGLE version number.
set(EAGLE_VERSION_MAJOR "2")
set(EAGLE_VERSION_MINOR "0")
set(EAGLE_VERSION_PATCH "3")
set(EAGLE_VERSION_DEVEL "")
set(EAGLE_VERSION "${EAGLE_VERSION_MAJOR}.${EAGLE_VERSION_MINOR}.${EAGLE_VERSION_PATCH}${EAGLE_VERSION_DEVEL}")
set(EAGLE_VERSION_FULL "${PROJECT_NAME}-${EAGLE_VERSION}")

if (CPACK_GENERATOR)
    message (STATUS "Configuring to produce the following package types: ${CPACK_GENERATOR}")
    SET(CPACK_PACKAGE_VENDOR "Illumina")
    SET(CPACK_PACKAGE_VERSION_MAJOR "${EAGLE_VERSION_MAJOR}")
    SET(CPACK_PACKAGE_VERSION_MINOR "${EAGLE_VERSION_MINOR}")
    SET(CPACK_PACKAGE_VERSION_PATCH "${EAGLE_VERSION_PATCH}${EAGLE_VERSION_DEVEL}")
    SET(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
    SET(CPACK_SET_DESTDIR ON)
    INCLUDE(CPack)
endif (CPACK_GENERATOR)

# Installation directories
if    (NOT EAGLE_PREFIX)
    set(EAGLE_PREFIX "/usr/local")
endif (NOT EAGLE_PREFIX)

if    (NOT CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "${EAGLE_PREFIX}")
endif (NOT CMAKE_INSTALL_PREFIX)

if (CMAKE_INSTALL_PREFIX)
    string(SUBSTRING "${CMAKE_INSTALL_PREFIX}" 0 1  slash)
    string(COMPARE NOTEQUAL "/" "${slash}" isslash)
    if    (isslash)
        install(CODE "
            string(LENGTH \"\$ENV{DESTDIR}\" ddlen)
            if (ddlen)
                message (FATAL_ERROR \"Installation prefix must begin with '/' if DESTDIR is set.\")
            endif (ddlen)
        ")
    endif (isslash)
endif (CMAKE_INSTALL_PREFIX)

if    (NOT EAGLE_EXEC_PREFIX)
    set(EAGLE_EXEC_PREFIX "." CACHE PATH "Install bin directory" FORCE)
endif (NOT EAGLE_EXEC_PREFIX)

if    (NOT EAGLE_BINDIR)
    set(EAGLE_BINDIR "${EAGLE_EXEC_PREFIX}/bin" 
        CACHE PATH "Install bin directory" FORCE)
endif (NOT EAGLE_BINDIR)

if    (NOT EAGLE_LIBDIR)
    set(EAGLE_LIBDIR "${EAGLE_EXEC_PREFIX}/lib/${EAGLE_VERSION_FULL}" 
        CACHE PATH "Install lib directory" FORCE)
endif (NOT EAGLE_LIBDIR)
set(EAGLE_PERL_LIBDIR "${EAGLE_LIBDIR}" CACHE PATH "Install perl lib directory" FORCE)

if    (NOT EAGLE_LIBEXECDIR)
    set(EAGLE_LIBEXECDIR "${EAGLE_EXEC_PREFIX}/libexec/${EAGLE_VERSION_FULL}" 
        CACHE PATH "Install libexec directory" FORCE)
endif (NOT EAGLE_LIBEXECDIR)

if    (NOT EAGLE_INCLUDEDIR)
    set(EAGLE_INCLUDEDIR "include/${EAGLE_VERSION_FULL}" 
        CACHE PATH "Install include directory" FORCE)
endif (NOT EAGLE_INCLUDEDIR)

if    (NOT EAGLE_DATADIR)
    set(EAGLE_DATADIR "share/${EAGLE_VERSION_FULL}" 
        CACHE PATH "Install data directory" FORCE)
endif (NOT EAGLE_DATADIR)

if    (NOT EAGLE_ETCDIR)
    set(EAGLE_ETCDIR "${EAGLE_EXEC_PREFIX}/etc/${EAGLE_VERSION_FULL}" 
        CACHE PATH "Install sysconfdir directory" FORCE)
endif (NOT EAGLE_ETCDIR)

if    (NOT EAGLE_DOCDIR)
    set(EAGLE_DOCDIR "doc/${EAGLE_VERSION_FULL}" 
        CACHE PATH "Install doc directory" FORCE)
endif (NOT EAGLE_DOCDIR)

if    (NOT EAGLE_MANDIR)
    set(EAGLE_MANDIR "man" 
        CACHE PATH "Install man directory" FORCE)
endif (NOT EAGLE_MANDIR)

message (STATUS "install        prefix: ${EAGLE_PREFIX}")
message (STATUS "install   exec prefix: ${EAGLE_EXEC_PREFIX}")
message (STATUS "install       bin dir: ${EAGLE_BINDIR}")
message (STATUS "install       lib dir: ${EAGLE_LIBDIR}")
message (STATUS "install   libexec dir: ${EAGLE_LIBEXECDIR}")
message (STATUS "install   include dir: ${EAGLE_INCLUDEDIR}")
message (STATUS "install      data dir: ${EAGLE_DATADIR}")
message (STATUS "install       doc dir: ${EAGLE_DOCDIR}")
message (STATUS "install       man dir: ${EAGLE_MANDIR}")

add_custom_target(EAGLE_CXX)

set(BOOST_REDIST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../redist")

# required boost libraries
set      (EAGLE_BOOST_VERSION 1.47.0)
set      (EAGLE_BOOST_COMPONENTS date_time filesystem iostreams program_options
                                regex serialization system thread)
set      (Boost_USE_MULTITHREAD ON)

set (EAGLE_GLOBALS_CMAKE "${CMAKE_SOURCE_DIR}/cmake/globals.cmake")
set (EAGLE_MACROS_CMAKE "${CMAKE_SOURCE_DIR}/cmake/macros.cmake")

include ("${CMAKE_SOURCE_DIR}/cmake/boost.cmake")

# required libzoo libraries
set      (LIBZOO_REDIST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../redist")
set      (EAGLE_LIBZOO_VERSION 1.0)
include ("${CMAKE_SOURCE_DIR}/cmake/libzoo.cmake")

add_subdirectory (cmake/preInstall)
add_subdirectory (bash)
add_subdirectory (c++)
add_subdirectory (perl)
#add_subdirectory (python)
add_subdirectory (makefiles)
#add_subdirectory (templates)
add_subdirectory ("${CMAKE_SOURCE_DIR}/../data" "${EAGLE_DATADIR}")
add_subdirectory ("${CMAKE_SOURCE_DIR}/../validation" "${CMAKE_BINARY_DIR}/validation")
add_subdirectory (cmake/postInstall)

if(BUILD_TESTING)
    set (CTEST_DATA_DIR "${CMAKE_SOURCE_DIR}/../data/tests")
    set (CTEST_SCRIPTS_DIR "${CMAKE_BINARY_DIR}/Testing/scripts")
    set (EAGLE_ALL_TESTEE applyVariants)
    include ("${CMAKE_SOURCE_DIR}/cmake/validation.cmake")
endif(BUILD_TESTING)
