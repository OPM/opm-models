#!/bin/bash

source `dirname $0`/build-ewoms.sh

declare -a upstreams
upstreams=(opm-parser
           opm-material
           opm-core
           opm-grid)

declare -A upstreamRev
upstreamRev[opm-parser]=master
upstreamRev[opm-material]=master
upstreamRev[opm-core]=master
upstreamRev[opm-grid]=master

ERT_REVISION=master
OPM_COMMON_REVISION=master

build_ewoms
test $? -eq 0 || exit 1

cp serial/build-ewoms/testoutput.xml .
