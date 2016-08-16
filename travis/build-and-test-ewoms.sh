#!/usr/bin/env bash
set -e

pushd . > /dev/null
ewoms/travis/build-ewoms.sh
cd ewoms/build
ctest --output-on-failure
popd > /dev/null
