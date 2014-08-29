#!/bin/bash
################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file installBoost.sh
##
## Script to install boost
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
VERSION=`echo ${EAGLE_BOOST_VERSION} | sed "s/\./_/g"`
SOURCE_TARBALL=${REDIST_DIR}/boost_${VERSION}.tar.bz2
TARBALL_COMPRESSION=j
SOURCE_DIR=${BUILD_DIR}/boost_${VERSION}

common_options $@

if [[ $CLEAN ]] ; then
    echo removing $SOURCE_DIR
    rm -rf $SOURCE_DIR ${INCLUDE_DIR}/boost ${LIB_DIR}/libboost_*.{a,so}
    exit 0
fi

common_create_source
cd ${SOURCE_DIR} \
    && ./bootstrap.sh ${BOOTSTRAP_OPTIONS} --prefix=${INSTALL_DIR} --with-libraries=`echo ${EAGLE_BOOST_COMPONENTS} | sed "s/;/,/g"` \
    && ./bjam -j$PARALLEL ${BJAM_OPTIONS} --libdir=${INSTALL_DIR}/lib --layout=system link=static threading=multi install

if [ $? != 0 ] ; then echo "$SCRIPT: build failed: Terminating..." >&2 ; exit 1 ; fi

#echo "Cleaning up ${SOURCE_DIR}"  >&2
#rm -rf ${SOURCE_DIR}

echo "boost-$VERSION installed successfully"  >&2
