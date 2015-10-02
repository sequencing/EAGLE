#!/bin/bash
set -o nounset

catcher()
{
    echo "Error detected"
    exit 1
}
trap catcher ERR
trap catcher SIGSEGV
