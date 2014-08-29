################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## file common.sh
##
## Definition of functions and variables common to all bootstrap scripts.
##
## author Mauricio Varea
##
################################################################################

function common_options () {
    TEMP=`getopt -n $SCRIPT -o fc -- "$@"`
    if [ $? != 0 ] ; then echo $SCRIPT: invalid option  >&2; echo "Terminating..." >&2 ; exit 2 ; fi
    eval set -- "$TEMP"
    FORCE=
    CLEAN=
    while true ; do
        case "$1" in
            -f) FORCE=true ; shift ;;
            -c) CLEAN=true ; shift ;;
            --)              shift ; break ;;
            *) echo "Internal error!" >&2; exit 2 ;;
        esac
    done
}

function common_create_source () {
    if [[ ! -e $SOURCE_TARBALL ]] ; then
        echo $SCRIPT: source tarball $SOURCE_TARBALL not found >&2
        exit 2
    fi  
    echo Decompressing $SOURCE_TARBALL >&2
    mkdir -p ${BUILD_DIR}
    tar -C ${BUILD_DIR} -${TARBALL_COMPRESSION}xf $SOURCE_TARBALL
    
    if [[ ! -d $SOURCE_DIR ]] ; then
        echo $SOURCE_DIR does not exist >&2
        exit 2
    fi
}

