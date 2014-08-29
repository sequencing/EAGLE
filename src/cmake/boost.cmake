################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for boost installation
##
## author Mauricio Varea
##
################################################################################

macro (resetFindBoost)
    unset (Boost_FOUND CACHE)
    unset (Boost_INCLUDE_DIRS CACHE)
    unset (Boost_INCLUDE_DIR CACHE)
    unset (Boost_LIBRARIES CACHE)
    unset (Boost_LIBRARY_DIRS CACHE)
    unset (Boost_VERSION CACHE)
    unset (Boost_LIB_VERSION CACHE)
    unset (Boost_MAJOR_VERSION CACHE)
    unset (Boost_MINOR_VERSION CACHE)
    unset (Boost_SUBMINOR_VERSION CACHE)
    unset (Boost_USE_STATIC_LIBS CACHE) 

    unset (ENV{BOOST_LIBRARYDIR})
    unset (Boost_USE_MULTITHREADED CACHE)

    foreach (COMPONENT ${EAGLE_BOOST_COMPONENTS})
        STRING(TOUPPER ${COMPONENT} UPPERCOMPONENT)
        unset (Boost_${UPPERCOMPONENT}_FOUND CACHE)
        unset (Boost_${UPPERCOMPONENT}_LIBRARY CACHE)
        unset (Boost_${UPPERCOMPONENT}_LIBRARY_RELEASE CACHE)
        unset (Boost_${UPPERCOMPONENT}_LIBRARY_DEBUG CACHE)
    endforeach (COMPONENT ${EAGLE_BOOST_COMPONENTS})
    

    unset (Boost_FOUND)
    unset (Boost_INCLUDE_DIRS)
    unset (Boost_INCLUDE_DIR)
    unset (Boost_LIBRARIES)
    unset (Boost_LIBRARY_DIRS)
    unset (Boost_VERSION)
    unset (Boost_LIB_VERSION)
    unset (Boost_MAJOR_VERSION)
    unset (Boost_MINOR_VERSION)
    unset (Boost_SUBMINOR_VERSION)
    unset (Boost_USE_STATIC_LIBS)

    foreach (COMPONENT ${EAGLE_BOOST_COMPONENTS})
        STRING(TOUPPER ${COMPONENT} UPPERCOMPONENT)
        unset (Boost_${UPPERCOMPONENT}_FOUND)
        unset (Boost_${UPPERCOMPONENT}_LIBRARY)
        unset (Boost_${UPPERCOMPONENT}_LIBRARY_RELEASE)
        unset (Boost_${UPPERCOMPONENT}_LIBRARY_DEBUG)
    endforeach (COMPONENT ${EAGLE_BOOST_COMPONENTS})

endmacro (resetFindBoost)

#   
# Not only finds boost but also sets the variables so that
# it is being used for include and linking
# Also makes sure pthread is available for boost
#
macro(eagle_find_boost boost_version boost_components)

    # pthread library required by boost
    eagle_find_library(PTHREAD "pthread.h" pthread)
    if    (HAVE_PTHREAD)
        set  (EAGLE_ADDITIONAL_LIB ${EAGLE_ADDITIONAL_LIB} pthread)
        message(STATUS "pthread supported")
    else  (HAVE_PTHREAD)
        message(STATUS "pthread headers: ${PTHREAD_INCLUDE_DIR}")
        message(STATUS "pthread library: ${PTHREAD_LIBRARY}")
        message(FATAL_ERROR "pthread library is required to build the EAGLE simulator")
    endif (HAVE_PTHREAD)

    find_package(Boost ${boost_version} REQUIRED ${boost_components})

    include_directories(BEFORE SYSTEM ${Boost_INCLUDE_DIRS})

    set      (HAVE_LIBBOOST_DATE_TIME       ${Boost_DATE_TIME_FOUND})
    set      (HAVE_LIBBOOST_FILESYSTEM      ${Boost_FILESYSTEM_FOUND})
    set      (HAVE_LIBBOOST_IOSTREAMS       ${Boost_IOSTREAMS_FOUND})
    set      (HAVE_LIBBOOST_PROGRAM_OPTIONS ${Boost_PROGRAM_OPTIONS_FOUND})
    set      (HAVE_LIBBOOST_PYTHON          ${Boost_PYTHON_FOUND})
    set      (HAVE_LIBBOOST_REGEX           ${Boost_REGEX_FOUND})
    set      (HAVE_LIBBOOST_SERIALIZATION   ${Boost_SERIALIZATION_FOUND})
    set      (HAVE_LIBBOOST_SYSTEM          ${Boost_SYSTEM_FOUND})
endmacro(eagle_find_boost)


if (EAGLE_FORCE_STATIC_LINK)
    set(Boost_USE_STATIC_LIBS ON) 
endif (EAGLE_FORCE_STATIC_LINK)

find_package(Boost ${EAGLE_BOOST_VERSION} COMPONENTS ${EAGLE_BOOST_COMPONENTS})
#
# If the right version of boost is not found, it will be built from the distribution
#
if (NOT Boost_FOUND)
    if (BOOSTROOT)
        message (STATUS "BOOSTROOT is set to ${BOOSTROOT} but boost ${EAGLE_BOOST_VERSION} was not found.")
        message (FATAL_ERROR "Unset BOOSTROOT or set it to the root location of boost ${EAGLE_BOOST_VERSION}.")
    endif(BOOSTROOT)
    if (BOOST_ROOT)
        message (STATUS "BOOST_ROOT is set to ${BOOST_ROOT} but boost ${EAGLE_BOOST_VERSION} was not found.")
        message (FATAL_ERROR "Unset BOOST_ROOT or set it to the root location of boost ${EAGLE_BOOST_VERSION}.")
    endif(BOOST_ROOT)

    # Try to find it in target installation location
    resetFindBoost()
    message(FATAL_ERROR "Boost ${EAGLE_BOOST_VERSION} not found. You need to set BOOST_ROOT to the root location of boost.")

endif (NOT Boost_FOUND)

