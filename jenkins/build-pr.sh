#!/bin/bash

source `dirname $0`/build-ewoms.sh

declare -a upstreams
upstreams=(opm-parser
           opm-output
           opm-material
           opm-core
           opm-grid)

declare -A upstreamRev
upstreamRev[opm-parser]=master
upstreamRev[opm-output]=master
upstreamRev[opm-material]=master
upstreamRev[opm-core]=master
upstreamRev[opm-grid]=master

ERT_REVISION=master
OPM_COMMON_REVISION=master

if grep -q "ert=" <<< $ghprbCommentBody
then
  ERT_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*ert=([0-9]+).*/\1/g'`/merge
fi

if grep -q "opm-common=" <<< $ghprbCommentBody
then
  OPM_COMMON_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*opm-common=([0-9]+).*/\1/g'`/merge
fi

for upstream in ${upstreams[*]}
do
  if grep -q "$upstream=" <<< $ghprbCommentBody
  then
    upstreamRev[$upstream]=pull/`echo $ghprbCommentBody | sed -r "s/.*$upstream=([0-9]+).*/\1/g"`/merge
  fi
done

echo "Building with ert=$ERT_REVISION opm-common=$OPM_COMMON_REVISION opm-parser=${upstreamRev[opm-parser]} opm-output=${upstreamRev[opm-output]} opm-material=${upstreamRev[opm-material]} opm-core=${upstreamRev[opm-core]} opm-grid=${upstreamRev[opm-grid]} ewoms=$sha1"

build_ewoms

test $? -eq 0 || exit 1

cp serial/build-ewoms/testoutput.xml .
