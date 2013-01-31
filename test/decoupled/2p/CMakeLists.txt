# set the CMake module include path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# eWoms unit-testing infrastructure
include(EwomsAddTest)

EwomsAddTest(test_impes)
EwomsAddTest(test_impesadaptive CONDITION ALUGRID_FOUND)
EwomsAddTest(test_transport)

# the MPFA tests are a bit sketchy: all tests are multiplexed through
# a single executable, and the name of this executable does not
# correspond to any test name, this means that we have to compile the
# executable ourselfs and jump though quite a few hoops to make the
# tests work with the eWoms testing infrastructure macro.
EwomsAddTest(test_mpfa2p ONLY_COMPILE)
EwomsAddTest(test_mpfao2p
  NO_COMPILE
  CONDITION ALUGRID_FOUND
  EXE_NAME test_mpfa2p
  DRIVER_ARGS --simulation test_mpfa2p -ModelType MPFAO)
EwomsAddTest(test_mpfal2p
  NO_COMPILE
  CONDITION ALUGRID_FOUND
  EXE_NAME test_mpfa2p
  DRIVER_ARGS --simulation test_mpfa2p -ModelType MPFAL)
EwomsAddTest(test_mpfal2padaptive
  NO_COMPILE
  CONDITION ALUGRID_FOUND
  EXE_NAME test_mpfa2p
  DRIVER_ARGS --simulation test_mpfa2p -ModelType MPFALAdaptive --end-time=2e4)
