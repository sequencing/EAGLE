#!/usr/bin/env bash
################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file check-source.sh
##
## Basic sanity checks on the source code before running the cppuint tests
##
## author Come Raczy
##
################################################################################

target=$1
shift
good=yes
for file in $* ; do
    check=`grep -nH CPPUNIT_TEST_SUITE_NAMED_REGISTRATION $file | grep -v registryName`
    if [[ $check ]] ; then
        if [[ $good ]] ; then
            good=
            echo >&2
            echo use of unchecked registry names: >&2
        fi
        echo "    "$check >&2
        echo >&2
    fi
done
if [[ $good ]] ; then
    echo checked > $target
else 
    exit 1
fi

