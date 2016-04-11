#!/bin/bash

source `dirname $0`/build-ewoms.sh

OPM_COMMON_REVISION=master
OPM_PARSER_REVISION=master
OPM_MATERIAL_REVISION=master
OPM_CORE_REVISION=master
OPM_GRID_REVISION=master
EWOMS_REVISION=master

build_ewoms
test $? -eq 0 || exit 1

cp serial/build-ewoms/testoutput.xml .
