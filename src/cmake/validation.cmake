################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file validation.cmake
##
## CMake configuration file for a validation framework
##
## authors: Mauricio Varea
##
################################################################################

macro(configure_tests dataDir scriptsDir suite bin)
    foreach(scenario ${suite})
        FILE(STRINGS "${dataDir}/${scenario}.txt" CTEST_LIST)
        message (STATUS "    Suite '${scenario}' includes the following integration tests:")
        message (STATUS "        ${CTEST_LIST}")
        foreach(TestId ${CTEST_LIST})
            add_test(NAME ${scenario}_${TestId} COMMAND ${scriptsDir}/test${scenario}.sh $<TARGET_FILE:${bin}> ${scenario} ${TestId})
        endforeach(TestId)
    endforeach(scenario)
endmacro(configure_tests)

foreach(eagleTool ${EAGLE_ALL_TESTEE})
    FILE(STRINGS "${CTEST_DATA_DIR}/${eagleTool}.txt" SUITE)
    message (STATUS "Validating tool '${eagleTool}' using the following test suites: '${SUITE}'")
    configure_tests("${CTEST_DATA_DIR}"  "${CTEST_SCRIPTS_DIR}"  "${SUITE}"  "${eagleTool}")
endforeach(eagleTool)

#  --test-action Start --test-action Build --test-action Test
add_custom_target(check  COMMAND ${CMAKE_CTEST_COMMAND} --test-model Experimental --test-action MemCheck --test-action Coverage )
add_custom_target(coverage-and-submit  COMMAND ${CMAKE_CTEST_COMMAND} --test-model Experimental --test-action Coverage --test-action Submit )
add_custom_target(memory-report COMMAND echo ''\; echo 'Memory Report Summary:' \;
                  find Testing -name \"test_?.log\" -exec awk '/ERROR SUMMARY/,FS=":" {printf \"%50s:%6d:%s\\n\", FILENAME, FNR, $$0}' {} + | sed 's/\(.*\)//g' | sed 's/ ERROR SUMMARY//g' )
add_custom_target(list-of-tests COMMAND ${CMAKE_CTEST_COMMAND} -N)

#add_dependencies( memory-report check )

