#!/bin/bash

function build_ewoms {
  # Build ERT
  pushd .
  mkdir -p $WORKSPACE/deps/ert
  cd $WORKSPACE/deps/ert
  git init .
  git remote add origin https://github.com/Ensembles/ert
  git fetch --depth 1 origin $ERT:branch_to_build
  test $? -eq 0 || exit 1
  git checkout branch_to_build
  popd

  pushd .
  mkdir -p serial/build-ert
  cd serial/build-ert
  cmake $WORKSPACE/deps/ert/devel -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install
  test $? -eq 0 || exit 1
  cmake --build . --target install
  test $? -eq 0 || exit 1
  popd

  # Build opm-common
  pushd .
  mkdir -p $WORKSPACE/deps/opm-common
  cd $WORKSPACE/deps/opm-common
  git init .
  git remote add origin https://github.com/OPM/opm-common
  git fetch origin $OPM_COMMON_REVISION:branch_to_build
  test $? -eq 0 || exit 1
  git checkout branch_to_build
  popd
  source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh

  pushd .
  mkdir serial/build-opm-common
  cd serial/build-opm-common
  build_module "-DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install -DCMAKE_PREFIX_PATH=$WORKSPACE/serial/install" 0 $WORKSPACE/deps/opm-common
  test $? -eq 0 || exit 1
  popd

  # Build opm-parser
  clone_and_build_module opm-parser "-DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install -DCMAKE_PREFIX_PATH=$WORKSPACE/serial/install" $OPM_PARSER_REVISION $WORKSPACE/serial
  test $? -eq 0 || exit 1

  # Build opm-material
  clone_and_build_module opm-material "-DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install -DCMAKE_PREFIX_PATH=$WORKSPACE/serial/install" $OPM_MATERIAL_REVISION $WORKSPACE/serial
  test $? -eq 0 || exit 1

  # Build opm-core
  clone_and_build_module opm-core "-DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install -DCMAKE_PREFIX_PATH=$WORKSPACE/serial/install" $OPM_CORE_REVISION $WORKSPACE/serial
  test $? -eq 0 || exit 1

  # Build opm-grid
  clone_and_build_module opm-grid "-DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install -DCMAKE_PREFIX_PATH=$WORKSPACE/serial/install" $OPM_GRID_REVISION $WORKSPACE/serial
  test $? -eq 0 || exit 1

  # Build ewoms
  git checkout $EWOMS_REVISION
  pushd .
  mkdir serial/build-ewoms
  cd serial/build-ewoms
  cmake -DCMAKE_PREFIX_PATH=$WORKSPACE/serial/install -DCMAKE_INSTALL_PREFIX=$WORKSPACE/serial/install -DUSE_QUADMATH=0 -DBUILD_TESTING=1 -DCMAKE_BUILD_TYPE=Release $WORKSPACE
  test $? -eq 0 || exit 1
  cmake --build .
  test $? -eq 0 || exit 1
  ctest -T Test --no-compress-output
  # Error code 8 appears to mean 'some tests skipped, but nothing failed'
  test $? -eq 8 || exit 1
  $WORKSPACE/deps/opm-common/jenkins/convert.py -x $WORKSPACE/deps/opm-common/jenkins/conv.xsl -t . > testoutput.xml
  popd
}
