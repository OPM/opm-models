# -*-cmake-*-
# - Try to find the DUNE localfunctions library
# Once done this will define:
#  DUNE_localfunctions_FOUND        - system has dune-localfunctions
#  DUNE_localfunctions_INCLUDE_DIR  - incude paths to use dune-localfunctions
#  DUNE_localfunctions_LIBRARIES    - Link these to use dune-localfunctions
INCLUDE(EwomsMacros)

EwomsSetup("DUNE_localfunctions" "dune-localfunctions" "DUNE")
EwomsParseDuneModuleInfo("dune-localfunctions" FILE_NAME "${DUNE_localfunctions_DIR}/dune.module")

EwomsFindIncludeDir("dune/localfunctions/lagrange.hh")

EwomsRequiredLibsFound()
EwomsIncludeDirsFound()
EwomsCheckFound()


