################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for the perl/lib directory
##
## author Mauricio Varea
##
################################################################################

include ("${EAGLE_GLOBALS_CMAKE}")

install(CODE "
    include (\"${EAGLE_MACROS_CMAKE}\")
    configure_files_recursively (\"${CMAKE_CURRENT_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}\" \"*.pm\")
    install_files_recursively (\"${CMAKE_CURRENT_BINARY_DIR}\" \"${EAGLE_ORIG_LIBDIR}\" \"*.pm\" \"\${EAGLE_LIBRARY_PERMISSIONS}\")
    ")
