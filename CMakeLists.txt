# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                              Name and version number  #
#       CMakeLists_files.cmake                   Path of source files     #
#       cmake/Modules/${project}-prereqs.cmake   Dependencies             #
#                                                                         #
###########################################################################

cmake_minimum_required (VERSION 2.8)

# additional search modules
find_path (OPM_MACROS_ROOT cmake/Modules/OpmInit.cmake
	PATHS ${opm-macros_ROOT}
	      ${CMAKE_CURRENT_SOURCE_DIR}
	DOC "Path to build system macros and scripts"
	NO_DEFAULT_PATH
	)
list (APPEND CMAKE_MODULE_PATH "${OPM_MACROS_ROOT}/cmake/Modules")

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include (${project}-prereqs)

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include (CMakeLists_files.cmake)

macro (config_hook)
	opm_need_version_of ("dune-common")
endmacro (config_hook)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (sources_hook)
endmacro (sources_hook)

macro (fortran_hook)
endmacro (fortran_hook)

macro (tests_hook)
endmacro (tests_hook)

# all setup common to the OPM library modules is done here
include (OpmLibMain)

file(GLOB grids "test/implicit/grids/*.dgf" "test/implicit/grids/*.art")
foreach(file ${grids})
  file(COPY ${file} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/grids)
endforeach()

include(EwomsAddTest)

# add targets for all tests of the models. we add the water-air test
# first because it take longest and so that we don't have to wait for
# them as long for parallel test runs
EwomsAddTest(waterair_pvs_ni SOURCES test/implicit/waterair_pvs_ni.cc)
EwomsAddTest(co2injection_immiscible SOURCES test/implicit/co2injection_immiscible.cc)
EwomsAddTest(co2injection_flash SOURCES test/implicit/co2injection_flash.cc)
EwomsAddTest(co2injection_ncp SOURCES test/implicit/co2injection_ncp.cc)
EwomsAddTest(co2injection_pvs SOURCES test/implicit/co2injection_pvs.cc)
EwomsAddTest(co2injection_immiscible_ni SOURCES test/implicit/co2injection_immiscible_ni.cc)
EwomsAddTest(co2injection_flash_ni SOURCES test/implicit/co2injection_flash_ni.cc)
EwomsAddTest(co2injection_ncp_ni SOURCES test/implicit/co2injection_ncp_ni.cc)
EwomsAddTest(co2injection_pvs_ni SOURCES test/implicit/co2injection_pvs_ni.cc)
EwomsAddTest(cuvette_pvs SOURCES test/implicit/cuvette_pvs.cc)
EwomsAddTest(diffusion_flash SOURCES test/implicit/diffusion_flash.cc)
EwomsAddTest(diffusion_ncp SOURCES test/implicit/diffusion_ncp.cc)
EwomsAddTest(diffusion_pvs SOURCES test/implicit/diffusion_pvs.cc)
EwomsAddTest(finger_immiscible SOURCES test/implicit/finger_immiscible.cc)
EwomsAddTest(fracture_discretefracture
  SOURCES test/implicit/fracture_discretefracture.cc
  CONDITION ALUGRID_FOUND
  DRIVER_ARGS --simulation fracture_discretefracture --end-time=1000)
EwomsAddTest(groundwater_immiscible SOURCES test/implicit/groundwater_immiscible.cc)
EwomsAddTest(infiltration_pvs SOURCES test/implicit/infiltration_pvs.cc)
EwomsAddTest(lens_immiscible SOURCES test/implicit/lens_immiscible.cc)
EwomsAddTest(lens_richards SOURCES test/implicit/lens_richards.cc)
EwomsAddTest(obstacle_immiscible SOURCES test/implicit/obstacle_immiscible.cc)
EwomsAddTest(obstacle_ncp SOURCES test/implicit/obstacle_ncp.cc)
EwomsAddTest(obstacle_pvs SOURCES test/implicit/obstacle_pvs.cc)
EwomsAddTest(outflow_pvs SOURCES test/implicit/outflow_pvs.cc)
EwomsAddTest(powerinjection_darcy SOURCES test/implicit/powerinjection_darcy.cc)
EwomsAddTest(powerinjection_forchheimer SOURCES test/implicit/powerinjection_forchheimer.cc)
EwomsAddTest(reservoir_blackoil SOURCES test/implicit/reservoir_blackoil.cc)
EwomsAddTest(test_navierstokes SOURCES test/implicit/test_navierstokes.cc CONDITION ALUGRID_FOUND AND SUPERLU_FOUND)
EwomsAddTest(test_quadrature SOURCES test/implicit/test_quadrature.cc DRIVER_ARGS --plain test_quadrature)
EwomsAddTest(test_stokes SOURCES test/implicit/test_stokes.cc CONDITION SUPERLU_FOUND)
EwomsAddTest(test_stokes2c SOURCES test/implicit/test_stokes2c.cc CONDITION SUPERLU_FOUND)
EwomsAddTest(test_stokesni SOURCES test/implicit/test_stokesni.cc CONDITION SUPERLU_FOUND)

EwomsAddTest(obstacle_immiscible_parallel
  SOURCES test/implicit/obstacle_immiscible.cc
  EXE_NAME obstacle_immiscible
  NO_COMPILE
  CONDITION MPI_FOUND
  DRIVER_ARGS --parallel-simulation obstacle_immiscible --end-time=1 --initial-time-step-size=1)

EwomsAddTest(obstacle_immiscible_parameters
  SOURCES test/implicit/obstacle_immiscible.cc
  EXE_NAME obstacle_immiscible
  NO_COMPILE
  DRIVER_ARGS --parameters obstacle_immiscible)

EwomsAddTest(obstacle_pvs_restart
  SOURCES test/implicit/obstacle_pvs.cc
  EXE_NAME obstacle_pvs
  NO_COMPILE
  DRIVER_ARGS --restart obstacle_pvs --pvs-verbosity=2 --end-time=12000)

# add the tutorial
EwomsAddTest(tutorial_coupled SOURCES tutorial/tutorial_coupled.cc)

# the tests for the common infrastructure
EwomsAddTest(test_propertysystem SOURCES test/common/propertysystem/test_propertysystem.cc DRIVER_ARGS --plain test_propertysystem)
