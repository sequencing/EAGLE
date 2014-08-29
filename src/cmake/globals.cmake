################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file globals.cmake
##
## CMake configuration file to identify the configuration of the system
##
## author Mauricio Varea
##
################################################################################

set(EAGLE_ORIG_ETCDIR      "${CMAKE_INSTALL_PREFIX}/${EAGLE_ETCDIR}")
set(EAGLE_ORIG_DATADIR     "${CMAKE_INSTALL_PREFIX}/${EAGLE_DATADIR}")
set(EAGLE_ORIG_BINDIR      "${CMAKE_INSTALL_PREFIX}/${EAGLE_BINDIR}")
set(EAGLE_ORIG_LIBDIR      "${CMAKE_INSTALL_PREFIX}/${EAGLE_LIBDIR}")
set(EAGLE_ORIG_LIBEXECDIR  "${CMAKE_INSTALL_PREFIX}/${EAGLE_LIBEXECDIR}")
set(EAGLE_ORIG_PERL_LIBDIR "${CMAKE_INSTALL_PREFIX}/${EAGLE_PERL_LIBDIR}")
set(EAGLE_ORIG_PYTHON_LIBDIR "${CMAKE_INSTALL_PREFIX}/${EAGLE_PYTHON_LIBDIR}")

install(CODE "

    # With package generator, the location where files are placed is not the location where they will be run.
    # _FULL_ variables are guaranteed valid only at runtime
    set (CPACK_GENERATOR \"${CPACK_GENERATOR}\")
    if (CPACK_GENERATOR)
        get_filename_component(EAGLE_FULL_ETCDIR       \"${EAGLE_ORIG_ETCDIR}\" ABSOLUTE)
        get_filename_component(EAGLE_FULL_DATADIR      \"${EAGLE_ORIG_DATADIR}\" ABSOLUTE)
        get_filename_component(EAGLE_FULL_BINDIR       \"${EAGLE_ORIG_BINDIR}\" ABSOLUTE)
        get_filename_component(EAGLE_FULL_LIBDIR       \"${EAGLE_ORIG_LIBDIR}\" ABSOLUTE)
        get_filename_component(EAGLE_FULL_LIBEXECDIR   \"${EAGLE_ORIG_LIBEXECDIR}\" ABSOLUTE)
        get_filename_component(EAGLE_FULL_PERL_LIBDIR  \"${EAGLE_ORIG_PERL_LIBDIR}\" ABSOLUTE)
        get_filename_component(EAGLE_FULL_PYTHON_LIBDIR  \"${EAGLE_ORIG_PYTHON_LIBDIR}\" ABSOLUTE)
    else (CPACK_GENERATOR)
        get_filename_component(EAGLE_FULL_ETCDIR       \"\$ENV{DESTDIR}${EAGLE_ORIG_ETCDIR}\" ABSOLUTE)
        get_filename_component(EAGLE_FULL_DATADIR      \"\$ENV{DESTDIR}${EAGLE_ORIG_DATADIR}\" ABSOLUTE)
        get_filename_component(EAGLE_FULL_BINDIR       \"\$ENV{DESTDIR}${EAGLE_ORIG_BINDIR}\" ABSOLUTE)
        get_filename_component(EAGLE_FULL_LIBDIR       \"\$ENV{DESTDIR}${EAGLE_ORIG_LIBDIR}\" ABSOLUTE)
        get_filename_component(EAGLE_FULL_LIBEXECDIR   \"\$ENV{DESTDIR}${EAGLE_ORIG_LIBEXECDIR}\" ABSOLUTE)
        get_filename_component(EAGLE_FULL_PERL_LIBDIR  \"\$ENV{DESTDIR}${EAGLE_ORIG_PERL_LIBDIR}\" ABSOLUTE)
        get_filename_component(EAGLE_FULL_PYTHON_LIBDIR  \"\$ENV{DESTDIR}${EAGLE_ORIG_PYTHON_LIBDIR}\" ABSOLUTE)
    endif (CPACK_GENERATOR)

    # _DEST_ variables always point to location where files are copied
    get_filename_component(EAGLE_DEST_ETCDIR       \"\$ENV{DESTDIR}${EAGLE_ORIG_ETCDIR}\" ABSOLUTE)
    get_filename_component(EAGLE_DEST_DATADIR      \"\$ENV{DESTDIR}${EAGLE_ORIG_DATADIR}\" ABSOLUTE)
    get_filename_component(EAGLE_DEST_BINDIR       \"\$ENV{DESTDIR}${EAGLE_ORIG_BINDIR}\" ABSOLUTE)
    get_filename_component(EAGLE_DEST_LIBDIR       \"\$ENV{DESTDIR}${EAGLE_ORIG_LIBDIR}\" ABSOLUTE)
    get_filename_component(EAGLE_DEST_LIBEXECDIR   \"\$ENV{DESTDIR}${EAGLE_ORIG_LIBEXECDIR}\" ABSOLUTE)
    get_filename_component(EAGLE_DEST_PERL_LIBDIR  \"\$ENV{DESTDIR}${EAGLE_ORIG_PERL_LIBDIR}\" ABSOLUTE)
    get_filename_component(EAGLE_DEST_PYTHON_LIBDIR  \"\$ENV{DESTDIR}${EAGLE_ORIG_PYTHON_LIBDIR}\" ABSOLUTE)

    set(EAGLE_VERSION_FULL \"${EAGLE_VERSION_FULL}\")
    set(EAGLE_VERSION \"${EAGLE_VERSION}\")
    set(EAGLE_DEBUG_MODE \"${EAGLE_DEBUG_MODE}\")
    set(EAGLE_EXECUTABLE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
    set(EAGLE_LIBRARY_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
    ")

