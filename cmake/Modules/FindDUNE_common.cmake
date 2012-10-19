# -*-cmake-*-
# - Try to find tje DUNE common library
# Once done this will define:
#  DUNE_common_FOUND        - system has dune-common
#  DUNE_common_INCLUDE_DIR  - incude paths to use dune-common
#  DUNE_common_LIBRARIES    - Link these to use dune-common
INCLUDE(EwomsMacros)

EwomsSetup("DUNE_common" "dune-common" "DUNE")

EwomsFindIncludeDir("dune/common/misc.hh")
EwomsFindLibrary("dunecommon")

EwomsRequiredLibsFound("dunecommon")
EwomsIncludeDirsFound()
EwomsCheckFound()


