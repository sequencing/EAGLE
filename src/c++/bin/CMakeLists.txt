################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file CMakeLists.txt
##
## Configuration file for the c++/bin subdirectory
##
## author Mauricio Varea
##
################################################################################

include(${EAGLE_CXX_EXECUTABLE_CMAKE})

file (GLOB EAGLE_PROGRAM_SOURCE_LIST [a-z][a-zA-Z0-9]*.cpp)

##
## Generic rule for all the other programs
##
#set(BAM_LIBRARY "${SAMTOOLS_DIR}/libbam.a")
foreach(EAGLE_PROGRAM_SOURCE ${EAGLE_PROGRAM_SOURCE_LIST})
    get_filename_component(EAGLE_PROGRAM ${EAGLE_PROGRAM_SOURCE} NAME_WE)
    add_executable        (${EAGLE_PROGRAM} ${EAGLE_PROGRAM_SOURCE})
    target_link_libraries (${EAGLE_PROGRAM} ${EAGLE_AVAILABLE_LIBRARIES}
                           ${BAM_LIBRARY} ${Boost_LIBRARIES} ${LIBXML2_LIBRARIES}  
                           ${EAGLE_ADDITIONAL_LIB} )
    install(TARGETS ${EAGLE_PROGRAM} RUNTIME DESTINATION ${EAGLE_BINDIR})
endforeach(EAGLE_PROGRAM_SOURCE)
