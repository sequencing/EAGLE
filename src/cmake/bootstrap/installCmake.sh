#!/bin/bash
################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file cwinstallCmakecommon.sh
##
## Installation script for cmake
##
## author Mauricio Varea
##
################################################################################

REDIST_DIR=$1
INSTALL_DIR=$2
if [[ $# -ge 3 ]] ; then PARALLEL=$3 ; else PARALLEL=1 ; fi

. `dirname "$0"`/common.sh

BUILD_DIR=${INSTALL_DIR}/build
BIN_DIR=${INSTALL_DIR}/bin
LIB_DIR=${INSTALL_DIR}/lib
INCLUDE_DIR=${INSTALL_DIR}/include

CMAKE_MAJOR=2
CMAKE_MINOR=8
CMAKE_PATCH=0
CMAKE_REQUIRED="$CMAKE_MAJOR.$CMAKE_MINOR.$CMAKE_PATCH"
TARBALL_VERSION="$CMAKE_MAJOR.$CMAKE_MINOR.4"
SCRIPT=`basename "$0"`
SOURCE_TARBALL=${REDIST_DIR}/cmake-$TARBALL_VERSION.tar.gz
TARBALL_COMPRESSION=z
SOURCE_DIR=${BUILD_DIR}/cmake-$TARBALL_VERSION
CMAKE_DIR=cmake-$CMAKE_MAJOR.$CMAKE_MINOR

common_options $@

check="^cmake version ([0-9]+)\.([0-9]+)\.([0-9]+)"
if [[ `cmake --version 2> /dev/null` =~ $check && ! $FORCE ]] ; then
    MAJOR=${BASH_REMATCH[1]}
    MINOR=${BASH_REMATCH[2]}
    PATCH=${BASH_REMATCH[3]}
    let VERSION_INT=1000000*MAJOR+1000*MINOR+PATCH
    let WANTED_VERSION_INT=1000000*CMAKE_MAJOR+1000*CMAKE_MINOR+CMAKE_PATCH
    if [[ "$VERSION_INT" -ge "$WANTED_VERSION_INT" ]] ; then
        echo "${BASH_REMATCH[0]} (>= $CMAKE_REQUIRED) detected" >&2
        exit 1
    fi
fi 

OLD_CMAKE_VERSION=`${BIN_DIR}/cmake --version 2> /dev/null`;
if [[ $OLD_CMAKE_VERSION == "cmake version $TARBALL_VERSION" && ! $FORCE ]] ; then
    echo cmake version \"$TARBALL_VERSION\" detected at ${BIN_DIR}/cmake >&2
    exit 0
elif [[ $OLD_CMAKE_VERSION != "" ]] ; then
    echo unable to install cmake version \"$TARBALL_VERSION\" in ${BIN_DIR} >&2 
    echo cmake version \"$OLD_CMAKE_VERSION\" is in the way. >&2
    echo Please use an empty location to build the product. >&2
    exit 2
fi 

echo "Cmake version >= 2.8 is required."
exit 3
