# -*-cmake-*-
# - Try to find tje DUNE istl library
# Once done this will define:
#  DUNE_istl_FOUND        - system has dune-istl
#  DUNE_istl_INCLUDE_DIR  - incude paths to use dune-istl
#  DUNE_istl_LIBRARIES    - Link these to use dune-istl
INCLUDE(EwomsMacros)

EwomsSetup("DUNE_istl" "dune-istl" "DUNE")
EwomsParseDuneModuleInfo("dune-istl" FILE_NAME "${DUNE_istl_DIR}/dune.module")

EwomsFindIncludeDir("dune/istl/io.hh")

EwomsRequiredLibsFound()
EwomsIncludeDirsFound()
EwomsCheckFound()


