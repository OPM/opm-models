#!/bin/bash

declare -a upstreams
upstreams=(opm-common
           ert
           opm-parser
           opm-output
           opm-material
           opm-core
           opm-grid)

declare -A upstreamRev
upstreamRev[opm-common]=master
upstreamRev[ert]=master
upstreamRev[opm-parser]=master
upstreamRev[opm-output]=master
upstreamRev[opm-material]=master
upstreamRev[opm-core]=master
upstreamRev[opm-grid]=master

if grep -q "opm-common=" <<< $ghprbCommentBody
then
  upstreamRev[opm-common]=pull/`echo $ghprbCommentBody | sed -r 's/.*opm-common=([0-9]+).*/\1/g'`/merge
fi

# Currently no downstream
declare -a downstreams
declare -A downstreamRev

# Clone opm-common
pushd .
mkdir -p $WORKSPACE/deps/opm-common
cd $WORKSPACE/deps/opm-common
git init .
git remote add origin https://github.com/OPM/opm-common
git fetch --depth 1 origin ${upstreamRev[opm-common]}:branch_to_build
test $? -eq 0 || exit 1
git checkout branch_to_build
popd

source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh

parseRevisions
printHeader ewoms

build_module_full ewoms
