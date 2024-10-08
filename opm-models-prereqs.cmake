# this avoids an annoying deprecation warning on DUNE 2.4 (which we
# are not interested in anyway)
set(DUNE_AVOID_CAPABILITIES_IS_PARALLEL_DEPRECATION_WARNING 1)

# defines that must be present in config.h for our headers
set (opm-models_CONFIG_VAR
  HAVE_VALGRIND
  HAVE_DUNE_COMMON
  HAVE_DUNE_GEOMETRY
  HAVE_DUNE_GRID
  HAVE_DUNE_LOCALFUNCTIONS
  HAVE_DUNE_ISTL
  HAVE_DUNE_ALUGRID
  HAVE_DUNE_FEM
  HAVE_ECL_INPUT
  HAVE_ECL_OUTPUT
  HAVE_OPM_GRID
  DUNE_AVOID_CAPABILITIES_IS_PARALLEL_DEPRECATION_WARNING
  HAVE_FLOATING_POINT_FROM_CHARS
  )

# dependencies
set (opm-models_DEPS
  # Need boost::test
  "Boost 1.44.0
    COMPONENTS unit_test_framework REQUIRED"
  # DUNE prerequisites
  "dune-common REQUIRED"
  "dune-geometry REQUIRED"
  "dune-grid REQUIRED"
  "dune-istl REQUIRED"
  "opm-common REQUIRED"
  "dune-localfunctions"
  "dune-alugrid"
  "dune-fem"
  "opm-grid"
  # valgrind client requests
  "Valgrind"
  # quadruple precision floating point calculations
  "QuadMath"
  )

find_package_deps(opm-models)
