################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file configure
##
## Configuration file for the share subdirectory
##
## author Mauricio Varea
##
################################################################################

#include ("${EAGLE_MACROS_CMAKE}")

#get_filename_component(EAGLE_CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message (STATUS "Verifying target directories access")

include ("${EAGLE_GLOBALS_CMAKE}")
install(
    CODE "set(EAGLE_DEST_DIRS \"\${EAGLE_DEST_DATADIR}\" 
                             \"\${EAGLE_DEST_BINDIR}\" \"\${EAGLE_DEST_LIBDIR}\" 
                             \"\${EAGLE_DEST_LIBEXECDIR}\")"
    SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/checkTargetPathWritable.cmake"
)
