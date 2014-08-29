################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for the perl subdirectory
##
## author Lilian Janin, Mauricio Varea
##
################################################################################

include ("${EAGLE_GLOBALS_CMAKE}")

install(CODE "
    include (\"${EAGLE_MACROS_CMAKE}\")
    configure_files_recursively (\"${CMAKE_CURRENT_SOURCE_DIR}/bin\" \"${CMAKE_CURRENT_BINARY_DIR}/bin\" \"[a-z][a-zA-Z0-9]*.pl\")
    install_files_recursively (\"${CMAKE_CURRENT_BINARY_DIR}/bin\" \"${EAGLE_ORIG_BINDIR}\" \"[a-z][a-zA-Z0-9]*.pl\" \"\${EAGLE_EXECUTABLE_PERMISSIONS}\")
    ")

##
## build all the libraries for the project
##
add_subdirectory (lib)

##
## build all the internal applications for the project
##
add_subdirectory (libexec)

## Temporarily:
#file (GLOB EAGLE_PROGRAM_LIST configure[a-zA-Z0-9]*.pl)
#install(FILES ${EAGLE_PROGRAM_LIST}
#        DESTINATION ${EAGLE_BINDIR}
#        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)

