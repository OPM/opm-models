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

# For eWoms, we set the BUILD_TESTING variable to off by default. The
# reason is that builing the tests takes quite some time and they are
# not very interesting for people who just want to use eWoms but not
# develop it. To enable testing, either pass --enable-testing to the
# configure script or call cmake with -DBUILD_TESTING=ON
set(MY_BUILD_TESTING "${BUILD_TESTING}")
set(BUILD_TESTING "ON" CACHE BOOL "Build the tests" FORCE)

# all setup common to the OPM library modules is done here
include (OpmLibMain)

file(GLOB datafiles "tests/data/*")
foreach(file ${datafiles})
  file(COPY ${file} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)
endforeach()

if (NOT MY_BUILD_TESTING)
	set(BUILD_TESTING "${MY_BUILD_TESTING}" CACHE BOOL "Build the tests" FORCE)
endif()

include(EwomsAddTest)

# add targets for all tests of the models. we add the water-air test
# first because it take longest and so that we don't have to wait for
# them as long for parallel test runs
EwomsAddTest(waterair_pvs_ni SOURCES tests/waterair_pvs_ni.cc)
EwomsAddTest(co2injection_immiscible_vcfv SOURCES tests/co2injection_immiscible_vcfv.cc)
EwomsAddTest(co2injection_flash_vcfv SOURCES tests/co2injection_flash_vcfv.cc)
EwomsAddTest(co2injection_ncp_vcfv SOURCES tests/co2injection_ncp_vcfv.cc)
EwomsAddTest(co2injection_pvs_vcfv SOURCES tests/co2injection_pvs_vcfv.cc)
EwomsAddTest(co2injection_immiscible_ni_vcfv SOURCES tests/co2injection_immiscible_ni_vcfv.cc)
EwomsAddTest(co2injection_flash_ni_vcfv SOURCES tests/co2injection_flash_ni_vcfv.cc)
EwomsAddTest(co2injection_ncp_ni_vcfv SOURCES tests/co2injection_ncp_ni_vcfv.cc)
EwomsAddTest(co2injection_pvs_ni_vcfv SOURCES tests/co2injection_pvs_ni_vcfv.cc)
EwomsAddTest(co2injection_immiscible_ecfv SOURCES tests/co2injection_immiscible_ecfv.cc)
EwomsAddTest(co2injection_flash_ecfv SOURCES tests/co2injection_flash_ecfv.cc)
EwomsAddTest(co2injection_ncp_ecfv SOURCES tests/co2injection_ncp_ecfv.cc)
EwomsAddTest(co2injection_pvs_ecfv SOURCES tests/co2injection_pvs_ecfv.cc)
EwomsAddTest(co2injection_immiscible_ni_ecfv SOURCES tests/co2injection_immiscible_ni_ecfv.cc)
EwomsAddTest(co2injection_flash_ni_ecfv SOURCES tests/co2injection_flash_ni_ecfv.cc)
EwomsAddTest(co2injection_ncp_ni_ecfv SOURCES tests/co2injection_ncp_ni_ecfv.cc)
EwomsAddTest(co2injection_pvs_ni_ecfv SOURCES tests/co2injection_pvs_ni_ecfv.cc)
EwomsAddTest(cuvette_pvs SOURCES tests/cuvette_pvs.cc)
EwomsAddTest(diffusion_flash SOURCES tests/diffusion_flash.cc)
EwomsAddTest(diffusion_ncp SOURCES tests/diffusion_ncp.cc)
EwomsAddTest(diffusion_pvs SOURCES tests/diffusion_pvs.cc)
EwomsAddTest(finger_immiscible SOURCES tests/finger_immiscible.cc CONDITION ALUGRID_FOUND)
EwomsAddTest(fracture_discretefracture
  SOURCES tests/fracture_discretefracture.cc
  CONDITION ALUGRID_FOUND
  DRIVER_ARGS --simulation fracture_discretefracture --end-time=1000)
EwomsAddTest(groundwater_immiscible SOURCES tests/groundwater_immiscible.cc)
EwomsAddTest(infiltration_pvs SOURCES tests/infiltration_pvs.cc)
EwomsAddTest(lens_immiscible SOURCES tests/lens_immiscible.cc)
EwomsAddTest(lens_richards SOURCES tests/lens_richards.cc)
EwomsAddTest(obstacle_immiscible SOURCES tests/obstacle_immiscible.cc)
EwomsAddTest(obstacle_ncp SOURCES tests/obstacle_ncp.cc)
EwomsAddTest(obstacle_pvs SOURCES tests/obstacle_pvs.cc)
EwomsAddTest(outflow_pvs SOURCES tests/outflow_pvs.cc)
EwomsAddTest(powerinjection_darcy SOURCES tests/powerinjection_darcy.cc)
EwomsAddTest(powerinjection_forchheimer SOURCES tests/powerinjection_forchheimer.cc)
EwomsAddTest(reservoir_blackoil SOURCES tests/reservoir_blackoil.cc)
EwomsAddTest(test_navierstokes SOURCES tests/test_navierstokes.cc CONDITION ALUGRID_FOUND AND SUPERLU_FOUND)
EwomsAddTest(test_quadrature SOURCES tests/test_quadrature.cc DRIVER_ARGS --plain test_quadrature)
EwomsAddTest(test_stokes SOURCES tests/test_stokes.cc CONDITION SUPERLU_FOUND)
EwomsAddTest(test_stokes2c SOURCES tests/test_stokes2c.cc CONDITION SUPERLU_FOUND)
EwomsAddTest(test_stokesni SOURCES tests/test_stokesni.cc CONDITION SUPERLU_FOUND)

# the test for the general-purpose Eclipse simulator 
EwomsAddTest(ecl_blackoil
  SOURCES tests/ecl_blackoil.cc
  EXE_NAME ecl_blackoil
  ONLY_COMPILE
  CONDITION HAVE_DUNE_CORNERPOINT_FOUND AND HAVE_ERT)

# test for the parallelization of the element centered finite volume
# discretization (using the non-isothermal NCP model and BiCGSTAB +
# ILU0)
EwomsAddTest(co2injection_ncp_ni_ecfv_parallel
  EXE_NAME co2injection_ncp_ni_ecfv
  NO_COMPILE
	PROCESSORS 4
  CONDITION MPI_FOUND
  DRIVER_ARGS --parallel-simulation=4 co2injection_ncp_ni_ecfv --end-time=10 --initial-time-step-size=1)

# test for the parallelization of the vertex centered finite volume
# discretization (using BiCGSTAB + ILU0)
EwomsAddTest(obstacle_immiscible_parallel
  EXE_NAME obstacle_immiscible
  NO_COMPILE
	PROCESSORS 4
  CONDITION MPI_FOUND
  DRIVER_ARGS --parallel-simulation=4 obstacle_immiscible --end-time=1 --initial-time-step-size=1)

# test for the parallel AMG linear solver using the vertex centered
# finite volume discretization
EwomsAddTest(lens_parallel
  EXE_NAME lens_immiscible
  NO_COMPILE
	PROCESSORS 4
  CONDITION MPI_FOUND
  DRIVER_ARGS --parallel-simulation=4 lens_immiscible --end-time=250 --initial-time-step-size=250)

EwomsAddTest(obstacle_immiscible_parameters
  EXE_NAME obstacle_immiscible
  NO_COMPILE
	DEPENDS obstacle_immiscible
  DRIVER_ARGS --parameters obstacle_immiscible)

EwomsAddTest(obstacle_pvs_restart
  EXE_NAME obstacle_pvs
  NO_COMPILE
	DEPENDS obstacle_pvs
  DRIVER_ARGS --restart obstacle_pvs --pvs-verbosity=2 --end-time=15000)

# add the tutorial
EwomsAddTest(tutorial1 SOURCES tutorial/tutorial1.cc)
