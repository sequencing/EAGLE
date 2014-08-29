################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file checkTargetPathWriteable.cmake
##
## Configuration file for the cppunit subfolders
##
## author Mauricio Varea
##
################################################################################

foreach (EAGLE_DEST_DIR ${EAGLE_DEST_DIRS})
    message (STATUS "Testing access to ${EAGLE_DEST_DIR}...")
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "make_directory" "${EAGLE_DEST_DIR}/.eagle" RESULT_VARIABLE TMP_RESULT )
    if (TMP_RESULT)
        message (STATUS "ERROR: Directory is not writeable: ${EAGLE_DEST_DIR}")
        message (STATUS "If you don't have administrator access to the "
                         "target installation location, please use --prefix "
                         "command-line option when configuring EAGLE. "
                         "Please use configure --help for all installer "
                         "command-line options details.")
        message (FATAL_ERROR "ERROR: EAGLE installation cannot continue")
    else (TMP_RESULT)
        execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "remove_directory" "${EAGLE_DEST_DIR}/.eagle" )
        message (STATUS "Directory is writeable: ${EAGLE_DEST_DIR}")
    endif (TMP_RESULT)
endforeach (EAGLE_DEST_DIR ${EAGLE_DEST_DIRS})
