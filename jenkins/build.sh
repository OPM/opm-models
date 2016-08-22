#!/bin/bash

source `dirname $0`/build-ewoms.sh

declare -a upstreams
upstreams=(ert
           opm-parser
           opm-output
           opm-material
           opm-core
           opm-grid)

declare -A upstreamRev
upstreamRev[ert]=master
upstreamRev[opm-parser]=master
upstreamRev[opm-output]=master
upstreamRev[opm-material]=master
upstreamRev[opm-core]=master
upstreamRev[opm-grid]=master

OPM_COMMON_REVISION=master

build_ewoms
test $? -eq 0 || exit 1

cp serial/build-ewoms/testoutput.xml .
