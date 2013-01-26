# set the CMake module include path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# eWoms unit-testing infrastructure
include(EwomsAddTest)

# we add the water-air test first because it take longest and so that
# we don't have to wait for them as long for parallel test runs
EwomsAddTest(waterair_pvs_ni)

# add targets for all tests of the implicit models
EwomsAddTest(co2injection_immiscible)
EwomsAddTest(co2injection_flash)
EwomsAddTest(co2injection_ncp)
EwomsAddTest(co2injection_pvs)
EwomsAddTest(co2injection_immiscible_ni)
EwomsAddTest(co2injection_flash_ni)
EwomsAddTest(co2injection_ncp_ni)
EwomsAddTest(co2injection_pvs_ni)
EwomsAddTest(cuvette_pvs)
EwomsAddTest(diffusion_flash)
EwomsAddTest(diffusion_ncp)
EwomsAddTest(diffusion_pvs)
EwomsAddTest(finger_immiscible)
EwomsAddTest(fracture_discretefracture CONDITION ALUGRID_FOUND)
EwomsAddTest(groundwater_immiscible)
EwomsAddTest(infiltration_pvs)
EwomsAddTest(lens_immiscible)
EwomsAddTest(lens_richards)
EwomsAddTest(obstacle_immiscible)
EwomsAddTest(obstacle_ncp)
EwomsAddTest(obstacle_pvs)
EwomsAddTest(outflow_pvs)
EwomsAddTest(powerinjection_darcy)
EwomsAddTest(powerinjection_forchheimer)
EwomsAddTest(reservoir_blackoil)
EwomsAddTest(test_navierstokes CONDITION ${ALUGRID_FOUND} AND ${SUPERLU_FOUND})
EwomsAddTest(test_quadrature DRIVER_ARGS --plain test_quadrature)
EwomsAddTest(test_stokes CONDITION SUPERLU_FOUND)
EwomsAddTest(test_stokes2c CONDITION SUPERLU_FOUND)
EwomsAddTest(test_stokesni CONDITION SUPERLU_FOUND)

EwomsAddTest(obstacle_immiscible_parallel
  EXE_NAME obstacle_immiscible
  NO_COMPILE
  CONDITION MPI_FOUND
  DRIVER_ARGS --parallel-simulation obstacle_immiscible --end-time=1 --initial-time-step-size=1)

EwomsAddTest(obstacle_immiscible_parameters
  EXE_NAME obstacle_immiscible
  NO_COMPILE
  DRIVER_ARGS --parameters obstacle_immiscible)

EwomsAddTest(obstacle_pvs_restart
  EXE_NAME obstacle_pvs
  NO_COMPILE
  DRIVER_ARGS --restart obstacle_pvs --pvs-verbosity=2)
