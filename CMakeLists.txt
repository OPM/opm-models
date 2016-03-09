# -*- mode: cmake; tab-width: 2; indent-tabs-mode: nil; truncate-lines: t
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

# add the current projects cmake module directory to the search
# path. This is not required anymore once support for federated builds
# is merged.
list(INSERT CMAKE_MODULE_PATH 0 "${PROJECT_SOURCE_DIR}/cmake/Modules")
set( OPM_COMMON_ROOT "" CACHE PATH "Root directory containing OPM related cmake modules")
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)

if(NOT OPM_COMMON_ROOT)
   find_package(opm-common QUIET)
endif()

if (opm-common_FOUND)
   include(OpmInit)
else()
   unset(opm-common_FOUND)

   if (NOT OPM_COMMON_ROOT AND SIBLING_SEARCH)
      set(OPM_COMMON_ROOT ${PROJECT_SOURCE_DIR}/../opm-common)
   endif()
   if (OPM_COMMON_ROOT)
      list( APPEND CMAKE_MODULE_PATH "${OPM_COMMON_ROOT}/cmake/Modules")
      include (OpmInit OPTIONAL RESULT_VARIABLE OPM_INIT)
      set( OPM_MACROS_ROOT ${OPM_COMMON_ROOT} )
   endif()

   if (NOT OPM_INIT)
      message( "" )
      message( " /---------------------------------------------------------------------------------\\")
      message( " |  Could not locate the opm build macros. The opm build macros                    |")
      message( " |  are in a separate repository - instructions to proceed:                        |")
      message( " |                                                                                 |")
      message( " |    1. Clone the repository: git clone git@github.com:OPM/opm-common.git         |")
      message( " |                                                                                 |")
      message( " |    2. Run cmake in the current project with -DOPM_COMMON_ROOT=<path>/opm-common |")
      message( " |                                                                                 |")
      message( " \\---------------------------------------------------------------------------------/")
      message( "" )
      message( FATAL_ERROR "Could not find OPM Macros")
   endif()

endif()

# Testing
include(OpmSatellites)

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include(${project}-prereqs)

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include(CMakeLists_files.cmake)

macro (config_hook)
  opm_need_version_of("dune-common")
  opm_need_version_of("dune-geometry")
  opm_need_version_of("dune-grid")
  opm_need_version_of("dune-localfunctions")
  opm_need_version_of("dune-alugrid")
  opm_need_version_of("dune-istl")
  opm_need_version_of("dune-fem")
endmacro (config_hook)

macro (files_hook)
endmacro (files_hook)

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

# all setup common to the OPM library modules
include (OpmLibMain)

# restore the value of the BUILD_TESTING variable to its original
# value. This is required because OpmLibMain thinks it is smarter
# appropriate and sets it...
if (NOT MY_BUILD_TESTING)
  set(BUILD_TESTING "${MY_BUILD_TESTING}" CACHE BOOL "Build the tests" FORCE)
endif()

# copy the data files (i.e., grid files, etc) needed to run the unit
# tests from the source to the build directory
file(GLOB datafiles "tests/data/*")
foreach(file ${datafiles})
  file(COPY ${file} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)
endforeach()

# include the macros required to build the tests and applications
# shipped as part of the eWoms distribution
include(EwomsAddApplication)

opm_set_test_driver("${PROJECT_SOURCE_DIR}/bin/runtest.sh" "--simulation")
opm_set_test_default_working_directory("${PROJECT_BINARY_DIR}")

set(CMAKE_PROJECT_NAME "${PROJECT_NAME}")

# the general-purpose ECL simulator ("ebos" == "&ecl &black-&oil &simulator)
EwomsAddApplication(ebos
                    SOURCES applications/ebos/ebos.cc
                    EXE_NAME ebos
                    CONDITION HAVE_OPM_GRID AND HAVE_ERT)

# the ART to DGF file format conversion utility
EwomsAddApplication(art2dgf
                    SOURCES applications/art2dgf.cc
                    EXE_NAME art2dgf)

opm_add_test(art2dgf
  NO_COMPILE
  DRIVER_ARGS --plain
  TEST_ARGS "data/fracture.art")

opm_add_test(test_blackoilfluidstate
  DRIVER_ARGS --plain)

# add targets for all tests of the models. we add the water-air test
# first because it take longest and so that we don't have to wait for
# them as long for parallel test runs
opm_add_test(waterair_pvs_ni
             TEST_ARGS --grid-global-refinements=1)

opm_add_test(lens_immiscible_vcfv
             TEST_ARGS --end-time=3000)

opm_add_test(lens_immiscible_ecfv
             TEST_ARGS --end-time=3000)

opm_add_test(finger_immiscible_ecfv
             CONDITION DUNE_ALUGRID_FOUND OR ALUGRID_FOUND)

opm_add_test(finger_immiscible_vcfv
             CONDITION DUNE_ALUGRID_FOUND OR ALUGRID_FOUND)

opm_add_test(finger_immiscible_ecfv_adaptive
             EXE_NAME finger_immiscible_ecfv
             CONDITION DUNE_ALUGRID_FOUND AND DUNE_FEM_FOUND
             NO_COMPILE
             DRIVER_ARGS --enable-grid-adaptation=true)

opm_add_test(test_navierstokes
             CONDITION (DUNE_ALUGRID_FOUND OR ALUGRID_FOUND) AND SUPERLU_FOUND)

opm_add_test(test_stokes2c
             CONDITION SUPERLU_FOUND)

opm_add_test(test_stokesni
             CONDITION SUPERLU_FOUND)

foreach(tapp co2injection_flash_ni_vcfv
             co2injection_flash_ni_ecfv
             co2injection_flash_vcfv
             co2injection_flash_ecfv
             co2injection_ncp_ni_vcfv
             co2injection_pvs_ni_vcfv
             co2injection_ncp_vcfv
             co2injection_pvs_vcfv
             co2injection_immiscible_ni_vcfv
             co2injection_immiscible_vcfv
             co2injection_immiscible_ecfv
             co2injection_ncp_ecfv
             co2injection_pvs_ecfv
             co2injection_immiscible_ni_ecfv
             co2injection_ncp_ni_ecfv
             co2injection_pvs_ni_ecfv
             powerinjection_forchheimer
             powerinjection_darcy
             cuvette_pvs
             infiltration_pvs
             lens_richards_vcfv
             lens_richards_ecfv
             obstacle_immiscible
             obstacle_ncp
             obstacle_pvs
             outflow_pvs
             diffusion_flash
             diffusion_ncp
             diffusion_pvs
             groundwater_immiscible)
  opm_add_test(${tapp})
endforeach()

opm_add_test(reservoir_blackoil_vcfv TEST_ARGS --end-time=8750000)
opm_add_test(reservoir_blackoil_ecfv TEST_ARGS --end-time=8750000)
opm_add_test(reservoir_ncp_vcfv TEST_ARGS --end-time=8750000)
opm_add_test(reservoir_ncp_ecfv TEST_ARGS --end-time=8750000)

opm_add_test(fracture_discretefracture
             CONDITION DUNE_ALUGRID_FOUND OR ALUGRID_FOUND
             TEST_ARGS --end-time=400)

opm_add_test(test_propertysystem
             DRIVER_ARGS --plain)

opm_add_test(test_quadrature
             DRIVER_ARGS --plain)

# test for the parallelization of the element centered finite volume
# discretization (using the non-isothermal NCP model and the parallel
# AMG linear solver)
opm_add_test(co2injection_ncp_ni_ecfv_parallel
             EXE_NAME co2injection_ncp_ni_ecfv
             NO_COMPILE
             PROCESSORS 4
             CONDITION MPI_FOUND
             DRIVER_ARGS --parallel-simulation=4)

# test for the parallelization of the vertex centered finite volume
# discretization (using BiCGSTAB + ILU0)
opm_add_test(obstacle_immiscible_parallel
             EXE_NAME obstacle_immiscible
             NO_COMPILE
             PROCESSORS 4
             CONDITION MPI_FOUND
             DRIVER_ARGS --parallel-simulation=4
             TEST_ARGS --end-time=1 --initial-time-step-size=1)

# test for the parallel AMG linear solver using the vertex centered
# finite volume discretization
opm_add_test(lens_immiscible_vcfv_parallel
             EXE_NAME lens_immiscible_vcfv
             NO_COMPILE
             PROCESSORS 4
             CONDITION MPI_FOUND
             DRIVER_ARGS --parallel-simulation=4
             TEST_ARGS --end-time=250 --initial-time-step-size=250)

opm_add_test(lens_immiscible_ecfv_parallel
             EXE_NAME lens_immiscible_ecfv
             NO_COMPILE
             PROCESSORS 4
             CONDITION MPI_FOUND
             DRIVER_ARGS --parallel-simulation=4
             TEST_ARGS --end-time=250 --initial-time-step-size=250)

opm_add_test(obstacle_immiscible_parameters
             EXE_NAME obstacle_immiscible
             NO_COMPILE
             DEPENDS obstacle_immiscible
             DRIVER_ARGS --parameters)

opm_add_test(obstacle_pvs_restart
             EXE_NAME obstacle_pvs
             NO_COMPILE
             DEPENDS obstacle_pvs
             DRIVER_ARGS --restart
             TEST_ARGS --pvs-verbosity=2 --end-time=15000)


# the test for the plain stokes model is completes quite
# quickly. thus, we put it to the end of the list...
opm_add_test(test_stokes
             CONDITION SUPERLU_FOUND)

opm_add_test(tutorial1
             SOURCES tutorial/tutorial1.cc)
