# -*-cmake-*-
# - Try to find the DUNE pdelab library
# Once done this will define:
#  DUNE_pdelab_FOUND        - system has dune-pdelab
#  DUNE_pdelab_INCLUDE_DIR  - incude paths to use dune-pdelab
#  DUNE_pdelab_LIBRARIES    - Link these to use dune-pdelab
INCLUDE(EwomsMacros)

EwomsSetup("DUNE_pdelab" "dune-pdelab" "DUNE")

EwomsFindIncludeDir("dune/pdelab/common/function.hh")

EwomsRequiredLibsFound()
EwomsIncludeDirsFound()
EwomsCheckFound()


