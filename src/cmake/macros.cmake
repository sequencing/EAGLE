################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file macros.cmake
##
## CMake configuration file for common installation macros
##
## authors: Mauricio Varea
##
################################################################################

macro(configure_files srcDir destDir pattern)
    file(GLOB templateFiles RELATIVE ${srcDir} ${srcDir}/${pattern}.in)
    foreach(templateFile ${templateFiles})
        string(REGEX REPLACE "\\.in" "" actualFile ${templateFile})
        message(STATUS "Configuring file ${srcDir}/${templateFile}") #...
        #message(STATUS "        ... into ${destDir}/${actualFile}")
        configure_file(${srcDir}/${templateFile} ${destDir}/${actualFile} @ONLY IMMEDIATE)
    endforeach(templateFile)
endmacro(configure_files)

macro(copy_dir srcDir destDir)
    file(GLOB_RECURSE templateFiles RELATIVE ${srcDir} ${srcDir}/*)
    foreach(templateFile ${templateFiles})
        message(STATUS "Copying file ${srcDir}/${templateFile}") #...
        #message(STATUS "    ... into ${destDir}/${templateFile}")
        configure_file(${srcDir}/${templateFile} ${destDir}/${templateFile} COPYONLY)
    endforeach(templateFile)
endmacro(copy_dir)

macro(install_files srcDir destDir pattern perm)
    file(GLOB templateFiles ${srcDir}/${pattern})
    file(INSTALL DESTINATION ${destDir} TYPE FILE
         FILES ${templateFiles} PERMISSIONS ${perm})
endmacro(install_files)

macro(configure_files_recursively srcDir destDir pattern)
    file(GLOB_RECURSE templateFiles RELATIVE ${srcDir} ${srcDir}/${pattern})
    foreach(templateFile ${templateFiles})
        message(STATUS "Configuring file ${srcDir}/${templateFile}")
        configure_file(${srcDir}/${templateFile} ${destDir}/${templateFile} @ONLY IMMEDIATE)
    endforeach(templateFile)
endmacro(configure_files_recursively)

macro(install_files_recursively srcDir destDir pattern perm)
    file(GLOB_RECURSE templateFiles RELATIVE ${srcDir} ${srcDir}/${pattern})
    foreach(templateFile ${templateFiles})
        get_filename_component(DIRNAME "${templateFile}" PATH)
        file(INSTALL DESTINATION ${destDir}/${DIRNAME} TYPE FILE
             FILES ${srcDir}/${templateFile} PERMISSIONS ${perm})
    endforeach(templateFile)
endmacro(install_files_recursively)

#   
# Macro to find libraries, with support for static-only search
#
macro(eagle_find_library name header library)
if    (NOT ${name}_INCLUDE_DIR)
    find_path(${name}_INCLUDE_DIR ${header} 
              HINTS ENV C_INCLUDE_PATH ENV CPATH ENV CPLUS_INCLUDE_PATH)
endif (NOT ${name}_INCLUDE_DIR)
if    (${name}_INCLUDE_DIR AND NOT ${name}_LIBRARY)
    find_library(${name}_LIBRARY 
                 NAMES "${EAGLE_LIBRARY_PREFIX}${library}${EAGLE_LIBRARY_SUFFIX}"
                 HINTS ENV LIBRARY_PATH)
endif (${name}_INCLUDE_DIR AND NOT ${name}_LIBRARY)
if(${name}_INCLUDE_DIR AND ${name}_LIBRARY)
    set (HAVE_${name} ${${name}_LIBRARY})
    message (STATUS "Found ${name}  header: ${${name}_INCLUDE_DIR}/${header}")
    message (STATUS "Found ${name} library: ${${name}_LIBRARY}")
endif(${name}_INCLUDE_DIR AND ${name}_LIBRARY)
endmacro(eagle_find_library)

