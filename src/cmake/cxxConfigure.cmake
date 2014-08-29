################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file cxxConfigure.cmake
##
## CMake configuration file for c++ executables
##
## author Mauricio Varea
##
################################################################################

INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(EAGLE_IS_BIG_ENDIAN)

INCLUDE(CheckFunctionExists)

find_path(HAVE_INTTYPES_H  inttypes.h)
find_path(HAVE_MEMORY_H    memory.h)
find_path(HAVE_STDINT_H    stdint.h)
find_path(HAVE_STDLIB_H    stdlib.h)
find_path(HAVE_STRING_H    string.h)
find_path(HAVE_STRINGS_H   strings.h)
find_path(HAVE_UNISTD_H    unistd.h)
find_path(HAVE_SYS_STAT_H  sys/stat.h)
find_path(HAVE_SYS_TYPES_H sys/types.h)

set (CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(floorf HAVE_FLOORF)
check_function_exists(round  HAVE_ROUND)
check_function_exists(roundf HAVE_ROUNDF)
check_function_exists(powf HAVE_POWF)
check_function_exists(erf HAVE_ERF)
check_function_exists(erf HAVE_ERFF)
check_function_exists(erfc HAVE_ERFC)
check_function_exists(erfc HAVE_ERFCF)

include ("${EAGLE_MACROS_CMAKE}")

# Support for static linking
# Note that this implies that all libraries must be found with the
# exact file name (libXXX.a or libXXX.so)
if    (EAGLE_FORCE_STATIC_LINK)
    message(STATUS "All libraries will be statically linked")
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-static")
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-static")
    # ensure that even if cmake decides to allow for dynamic libs resolution,
    # this gets overriden into static...
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS ${CMAKE_EXE_LINK_STATIC_CXX_FLAGS})
    set(EAGLE_LIBRARY_PREFIX ${CMAKE_STATIC_LIBRARY_PREFIX})
    set(EAGLE_LIBRARY_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
else  (EAGLE_FORCE_STATIC_LINK)
    set(EAGLE_LIBRARY_PREFIX "")
    set(EAGLE_LIBRARY_SUFFIX "")
endif (EAGLE_FORCE_STATIC_LINK)

if    (EAGLE_DEBUG_MODE)
    message(STATUS "Additional debugging mode:   ON")
else  (EAGLE_DEBUG_MODE)
    message(STATUS "Additional debugging mode:   OFF")
endif (EAGLE_DEBUG_MODE)

if    (EAGLE_SILENT_MODE)
    message(STATUS "Silent mode:                 ON")
else  (EAGLE_SILENT_MODE)
    message(STATUS "Silent mode:                 OFF")
endif (EAGLE_SILENT_MODE)

# optional support for gzip compression
eagle_find_library(ZLIB zlib.h z)
if    (HAVE_ZLIB)
    set  (EAGLE_ADDITIONAL_LIB ${EAGLE_ADDITIONAL_LIB} z)
    message(STATUS "gzip compression supported")
else  (HAVE_ZLIB)
    message(FATAL_ERROR "No support for gzip compression")
endif (HAVE_ZLIB)

# optional support for bzip2 compression
eagle_find_library(BZIP2 bzlib.h bz2)
if    (HAVE_BZIP2)
    set(HAVE_BZLIB HAVE_BZIP2)
    set(EAGLE_ADDITIONAL_LIB ${EAGLE_ADDITIONAL_LIB} bz2)
    message(STATUS "bzip2 compression supported")
else  (HAVE_BZIP2)
    message(FATAL_ERROR "No support for bzip2 compression")
endif (HAVE_BZIP2)

# optional support for libzoo
eagle_find_library(LIBZOO libzoo/cli/Common.hh zoo)
if    (HAVE_LIBZOO)
    include_directories(${LIBZOO_INCLUDE_DIR})
    link_libraries(${LIBZOO_LIBRARY})
    message(STATUS "libzoo supported")
else  (HAVE_LIBZOO)
    message(STATUS "No support for libzoo")
endif (HAVE_LIBZOO)

eagle_find_boost(${EAGLE_BOOST_VERSION} "${EAGLE_BOOST_COMPONENTS}")

eagle_find_library(CPGPLOT cpgplot.h cpgplot)
eagle_find_library(PGPLOT cpgplot.h pgplot)
eagle_find_library(X11 X.h X11)

include (FindLibXml2)
if    (NOT LIBXML2_FOUND)
    message (FATAL_ERROR "libxml2 was not found")
endif (NOT LIBXML2_FOUND)
include_directories(BEFORE SYSTEM ${LIBXML2_INCLUDE_DIR})

eagle_find_library(CPPUNIT "cppunit/Test.h" cppunit${CPPUNIT_DEBUG})

set (CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} -fopenmp -Wall -Wextra -Wunused -Wno-long-long -Wsign-compare -Wpointer-arith" CACHE STRING "g++ flags" FORCE)
set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "g++ flags" FORCE)
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "g++ flags" )
set (CMAKE_CXX_FLAGS_CDASH "-O0 -g -fprofile-arcs -ftest-coverage" CACHE STRING "g++ flags" )
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "g++ flags" FORCE)
set (CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG" CACHE STRING "g++ flags" FORCE)

# Force static linking
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE version)
    string(STRIP ${version} version)
    if    (version MATCHES "^4\\.0\\.2$")
        message (FATAL_ERROR "Unsupported GNU compiler version: "
                             "${version}: gcc bug \#28088 makes it impossible to "
                             "compile the EAGLE with gcc 4.0.2. "
                             "Please use gcc 4.1.2 instead.")
    endif (version MATCHES "^4\\.0\\.2$")
    
    string(REGEX REPLACE "^([0-9])\\.[0-9]\\.[0-9]" "\\1" major_version ${version})
    string(REGEX REPLACE "^[0-9]\\.([0-9])\\.[0-9]" "\\1" minor_version ${version})
    string(REGEX REPLACE "^[0-9]\\.[0-9]\\.([0-9])" "\\1" patch_version ${version})
    if    (major_version LESS 3 OR major_version EQUAL 3 AND minor_version LESS 4)
        message (FATAL_ERROR "Unsupported GNU C++ compiler: g++ version ${version}: "
                             "only g++ versions >= 3.4.0 are supported")
    endif (major_version LESS 3 OR major_version EQUAL 3 AND minor_version LESS 4)

    set("${CMAKE_CXX_COMPILER_ID}${major_version}" true)
    set("${CMAKE_CXX_COMPILER_ID}${major_version}${minor_version}" true)
    set("${CMAKE_CXX_COMPILER_ID}${major_version}${minor_version}${patch_version}" true)
    message (STATUS "using compiler: gcc version ${version}")

    if    (major_version EQUAL 3)
        message (STATUS "WARNING: there are known compatibility issues with "
                        "gcc versions 3.x. Please use gcc 4.x (except for "
                        "4.0.2) in case you experience linker failures.")
    endif (major_version EQUAL 3)

endif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

##
## Suppress spurious warnings in less recent compilers
##
if    (NOT GNU42)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter ")
endif (NOT GNU42)

if    (GNU412 OR GNU42 OR GNU43)
    ## Before 4.1.2, pedantic breaks on boost lambda expressions
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic ")
endif (GNU412 OR GNU42 OR GNU43)

if (CMAKE_SYSTEM_PROCESSOR MATCHES "^i[67]86$")
    ##
    ## Use scalar floating point instructions from the SSE instruction set.
    ## Note: Pentium3 SSE supports only single precision arithmetics
    ##
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -mfpmath=sse")
endif (CMAKE_SYSTEM_PROCESSOR MATCHES "^i[67]86$")
if (CMAKE_SYSTEM_PROCESSOR MATCHES "^i[345]86$")
    ##
    ## Prevent using 80bits registers (more consistent rounding)
    ##
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffloat-store")
endif (CMAKE_SYSTEM_PROCESSOR MATCHES "^i[345]86$")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lib/common/config.h.in ${EAGLE_CXX_CONFIG_H_DIR}/config.h)
