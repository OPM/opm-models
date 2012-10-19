# -*-cmake-*-
# - Try to find the Alberta grid manager
# Once done this will define:
#  Alberta_FOUND        - system has dune-grid
#  Alberta_INCLUDE_DIR  - incude paths to use dune-grid
#  Alberta_LIBRARIES    - Link these to use dune-grid
Include(EwomsMacros)

EwomsSetup("Alberta" "Alberta" "Alberta")

#EwomsAddPathSuffixes("${MyIncludeSuffixes}" "")

EwomsFindIncludeDir("alberta.h")
EwomsFindLibrary("ALBERTA22_0")
EwomsFindLibrary("ALBERTA22_1")
EwomsFindLibrary("alberta_util")

EwomsRequiredLibsFound("ALBERTA22_0" "ALBERTA22_1" "alberta_util")
EwomsIncludeDirsFound()
EwomsCheckFound()
