################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for libzoo installation
##
## author Lilian Janin
##
################################################################################


if (EAGLE_FORCE_STATIC_LINK)
    set(Libzoo_USE_STATIC_LIBS ON) 
endif (EAGLE_FORCE_STATIC_LINK)

#find_package(Libzoo ${EAGLE_LIBZOO_VERSION} COMPONENTS ${EAGLE_LIBZOO_COMPONENTS})
#
# If the right version of libzoo is not found, it will be built from the distribution
#
if (NOT Libzoo_FOUND)
    message(STATUS "Libzoo ${EAGLE_LIBZOO_VERSION} not found. Libzoo will be built from the distribution...")

    set(ENV{EAGLE_LIBZOO_VERSION} "${EAGLE_LIBZOO_VERSION}")
    if (NOT CMAKE_PARALLEL)
        set (CMAKE_PARALLEL "1")
    endif (NOT CMAKE_PARALLEL)
    execute_process(COMMAND "/bin/bash"
"${CMAKE_SOURCE_DIR}/cmake/bootstrap/installLibzoo.sh" "${LIBZOO_REDIST_DIR}"
"${CMAKE_CURRENT_BINARY_DIR}/bootstrap" "${CMAKE_PARALLEL}"  RESULT_VARIABLE TMP_RESULT )

    if (NOT TMP_RESULT)
        message(STATUS "Successfuly built libzoo ${EAGLE_LIBZOO_VERSION} from the distribution package...")
    else (NOT TMP_RESULT)
        message (FATAL_ERROR "Failed to build libzoo ${EAGLE_LIBZOO_VERSION}")
    endif (NOT TMP_RESULT)

    #set (LIBZOO_ROOT "${CMAKE_CURRENT_BINARY_DIR}/bootstrap")
    set (LIBZOO_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/bootstrap/build/beetl-1.0/src")
    set (LIBZOO_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/bootstrap/build/beetl-1.0/src/libzoo.a")

    #force static linking with redistributed libzoo.
    set (Libzoo_USE_STATIC_LIBS ON)

endif (NOT Libzoo_FOUND)

