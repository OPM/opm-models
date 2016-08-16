#!/usr/bin/env bash
set -e

pushd . > /dev/null
cd ewoms
mkdir build
cd build
cmake -DADD_DISABLED_CTESTS=OFF -DUSE_QUADMATH=OFF -DBUILD_TESTING=ON  ../
make
popd > /dev/null
