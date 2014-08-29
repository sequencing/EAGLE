#!/bin/bash
################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file installLibzoo.sh
##
## Script to install libzoo
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

SCRIPT=`basename "$0"`
VERSION=${EAGLE_LIBZOO_VERSION}
SOURCE_TARBALL=${REDIST_DIR}/beetl-${VERSION}.tar.bz2
TARBALL_COMPRESSION=j
SOURCE_DIR=${BUILD_DIR}/beetl-${VERSION}

common_options $@

if [[ $CLEAN ]] ; then
    echo removing $SOURCE_DIR
    rm -rf $SOURCE_DIR ${INCLUDE_DIR}/libzoo ${LIB_DIR}/liblibzoo_*.{a,so}
    exit 0
fi

common_create_source
cd ${SOURCE_DIR} \
    && ./configure --prefix=${INSTALL_DIR} \
    && make -j$PARALLEL -C src libzoo.a \
    && mkdir -p ${LIB_DIR} \
    && cp src/libzoo.a ${LIB_DIR}

if [ $? != 0 ] ; then echo "$SCRIPT: build failed: Terminating..." >&2 ; exit 1 ; fi

#echo "Cleaning up ${SOURCE_DIR}"  >&2
#rm -rf ${SOURCE_DIR}

echo "libzoo-$VERSION installed successfully"  >&2
