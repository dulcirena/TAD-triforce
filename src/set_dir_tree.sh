#!/usr/bin/bash

# Arg1: working directory
# Arg2: absolute path to the rawdata directory

wdir=$1
mkdir ${1}/out
mkdir ${1}/data
ln -s ${2}/* ${1}/data/.

