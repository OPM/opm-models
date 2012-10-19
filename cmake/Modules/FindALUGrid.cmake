# -*-cmake-*-
# - Try to find the UG grid manager
# Once done this will define:
#  ALUGrid_FOUND        - system has dune-grid
#  UG_INCLUDE_DIR  - incude paths to use dune-grid
#  UG_LIBRARIES    - Link these to use dune-grid
Include(EwomsMacros)

EwomsSetup("ALUGrid" "ALUGrid" "ALUGrid")

set(MyIncludeSuffixes 
    "include/serial"
    "include/parallel"
    "include/duneinterface")

EwomsAddPathSuffixes("${MyIncludeSuffixes}" "")

EwomsFindIncludeDir("alugrid_2d.h")
EwomsFindExtraIncludeDir("ALU_SERIAL" "serialize.h")
EwomsFindExtraIncludeDir("ALU_PARALLEL" "gitter_pll_impl.h")
EwomsFindExtraIncludeDir("ALU_DUNE" "gitter_dune_impl.h")
EwomsFindLibrary("alugrid")

EwomsRequiredLibsFound("alugrid")
EwomsIncludeDirsFound()
EwomsCheckFound()
