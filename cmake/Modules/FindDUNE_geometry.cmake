# -*-cmake-*-
# - Try to find the DUNE geometry library
# Once done this will define:
#  DUNE_geometry_FOUND        - system has dune-geometry
#  DUNE_geometry_INCLUDE_DIR  - incude paths to use dune-geometry
#  DUNE_geometry_LIBRARIES    - Link these to use dune-geometry
INCLUDE(EwomsMacros)

EwomsSetup("DUNE_geometry" "dune-geometry" "DUNE")

EwomsFindIncludeDir("dune/geometry/type.hh")
EwomsFindLibrary("dunegeometry")

EwomsRequiredLibsFound("dunegeometry")
EwomsIncludeDirsFound()
EwomsCheckFound()


