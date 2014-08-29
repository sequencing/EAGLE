################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for the test directory
##
## author Mauricio Varea
##
################################################################################

include ("${EAGLE_GLOBALS_CMAKE}")


# Install "test*"
file (GLOB EAGLE_TEST_LIST RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" TestConfig/test*.mk)

install(CODE "
    configure_file(\"${CMAKE_CURRENT_SOURCE_DIR}/Makefile.in\" \"${CMAKE_CURRENT_BINARY_DIR}/Makefile.top\" @ONLY)
    ")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Makefile.top" DESTINATION "${EAGLE_ORIG_DATADIR}/validation/" RENAME "Makefile")

foreach(EAGLE_TEST_CFG ${EAGLE_TEST_LIST})
    get_filename_component(EAGLE_TEST ${EAGLE_TEST_CFG} NAME_WE)
    message(STATUS "Configuring system-test: ${EAGLE_TEST}")
    install(CODE "
        set (EAGLE_TEST ${EAGLE_TEST})
        configure_file(\"${CMAKE_CURRENT_SOURCE_DIR}/MakefileForTests.in\" \"${CMAKE_CURRENT_BINARY_DIR}/Makefile.${EAGLE_TEST}\" @ONLY)
        ")
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Makefile.${EAGLE_TEST}" DESTINATION "${EAGLE_ORIG_DATADIR}/validation/${EAGLE_TEST}" RENAME "Makefile")
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/${EAGLE_TEST_CFG}" DESTINATION "${EAGLE_ORIG_DATADIR}/validation/${EAGLE_TEST}" RENAME "config.mk")
endforeach(EAGLE_TEST_CFG)


