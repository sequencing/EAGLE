################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for the data subdirectory
##
## author Mauricio Varea
##
################################################################################

include ("${EAGLE_GLOBALS_CMAKE}")
include ("${EAGLE_MACROS_CMAKE}")

# TODO: these ones could be generated from aggregating the 'acgt' dir
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/Genomes/acgt.fa" "${CMAKE_BINARY_DIR}/Testing/data/Genomes/acgt.fa" COPYONLY)
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/Genomes/acgt.fa.fai" "${CMAKE_BINARY_DIR}/Testing/data/Genomes/acgt.fa.fai" COPYONLY)
#
copy_dir( "${CMAKE_CURRENT_SOURCE_DIR}/Genomes/acgt" "${CMAKE_BINARY_DIR}/Testing/data/Genomes/acgt")
copy_dir( "${CMAKE_CURRENT_SOURCE_DIR}/tests" "${CMAKE_BINARY_DIR}/Testing/data")

install( FILES "${CMAKE_CURRENT_SOURCE_DIR}/Genomes/E_coli.fa" DESTINATION "${EAGLE_ORIG_DATADIR}/Genomes")
#install( FILES "${CMAKE_CURRENT_SOURCE_DIR}/RunInfo" DESTINATION "${EAGLE_ORIG_DATADIR}/")

install(CODE "
    include (\"${EAGLE_MACROS_CMAKE}\")
    install_files (\"${CMAKE_CURRENT_SOURCE_DIR}/QualityTables\" \"${EAGLE_ORIG_DATADIR}/QualityTables\" \"*\" \"\${EAGLE_LIBRARY_PERMISSIONS}\")
    install_files (\"${CMAKE_CURRENT_SOURCE_DIR}/MismatchTables\" \"${EAGLE_ORIG_DATADIR}/MismatchTables\" \"*\" \"\${EAGLE_LIBRARY_PERMISSIONS}\")
    install_files (\"${CMAKE_CURRENT_SOURCE_DIR}/HomopolymerQualityDropTables\" \"${EAGLE_ORIG_DATADIR}/HomopolymerQualityDropTables\" \"*\" \"\${EAGLE_LIBRARY_PERMISSIONS}\")
    install_files (\"${CMAKE_CURRENT_SOURCE_DIR}/MotifQualityDropTables\" \"${EAGLE_ORIG_DATADIR}/MotifQualityDropTables\" \"*\" \"\${EAGLE_LIBRARY_PERMISSIONS}\")
    install_files (\"${CMAKE_CURRENT_SOURCE_DIR}/TemplateLengthTables\" \"${EAGLE_ORIG_DATADIR}/TemplateLengthTables\" \"*\" \"\${EAGLE_LIBRARY_PERMISSIONS}\")
    install_files (\"${CMAKE_CURRENT_SOURCE_DIR}/GcCoverageFitTables\" \"${EAGLE_ORIG_DATADIR}/GcCoverageFitTables\" \"*\" \"\${EAGLE_LIBRARY_PERMISSIONS}\")
    install_files (\"${CMAKE_CURRENT_SOURCE_DIR}/RunInfo\" \"${EAGLE_ORIG_DATADIR}/RunInfo\" \"*\" \"\${EAGLE_LIBRARY_PERMISSIONS}\")
    install_files (\"${CMAKE_CURRENT_SOURCE_DIR}/Variants\" \"${EAGLE_ORIG_DATADIR}/Variants\" \"*\" \"\${EAGLE_LIBRARY_PERMISSIONS}\")
    execute_process(COMMAND rm -f \"${EAGLE_ORIG_DATADIR}/../EAGLE\"   )
    execute_process(COMMAND rm -f \"${EAGLE_ORIG_LIBDIR}/../EAGLE\"    )
    execute_process(COMMAND rm -f \"${EAGLE_ORIG_LIBEXECDIR}/../EAGLE\")
    execute_process(COMMAND ln -s \"${EAGLE_VERSION_FULL}\" \"${EAGLE_ORIG_DATADIR}/../EAGLE\"   )
    execute_process(COMMAND ln -s \"${EAGLE_VERSION_FULL}\" \"${EAGLE_ORIG_LIBDIR}/../EAGLE\"    )
    execute_process(COMMAND ln -s \"${EAGLE_VERSION_FULL}\" \"${EAGLE_ORIG_LIBEXECDIR}/../EAGLE\")
    ")
