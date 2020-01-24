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

# Mandatory call to project
project(opm-models C CXX)

cmake_minimum_required (VERSION 2.8)

# add the current projects cmake module directory to the search
# path. This is not required anymore once support for federated builds
# is merged.
list(INSERT CMAKE_MODULE_PATH 0 "${PROJECT_SOURCE_DIR}/cmake/Modules")
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)

if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
  else()
    string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
      # We are using build directories named <prefix><module-name><postfix>
      set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
    elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
      # All modules are in a common build dir
      set(opm-common_DIR "${_parent_full_dir}/opm-common")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)

include(OpmInit)

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
add_custom_target(opm-models_prepare)

# the ART to DGF file format conversion utility

if (BUILD_EXAMPLES)
  EwomsAddApplication(art2dgf
    SOURCES art2dgf/art2dgf.cc
    EXE_NAME art2dgf)

  add_dependencies(test-suite art2dgf)

  opm_add_test(art2dgf
    NO_COMPILE
    DRIVER_ARGS --plain
    TEST_ARGS "data/fracture-raw.art")
endif()

# add targets for all tests of the models. we add the water-air test
# first because it take longest and so that we don't have to wait for
# them as long for parallel test runs
opm_add_test(waterair_pvs_ni
             TEST_ARGS --grid-global-refinements=1)

opm_add_test(lens_immiscible_vcfv_ad
             TEST_ARGS --end-time=3000)

opm_add_test(lens_immiscible_vcfv_fd
             TEST_ARGS --end-time=3000)

opm_add_test(lens_immiscible_ecfv_ad
             TEST_ARGS --end-time=3000)

opm_add_test(lens_immiscible_ecfv_ad_23
             TEST_ARGS --end-time=3000)

# this test is identical to the simulation of the lens problem that
# uses the element centered finite volume discretization in
# conjunction with automatic differentiation
# (lens_immiscible_ecfv_ad). The only difference is that it uses
# multiple compile units in order to ensure that eWoms code can be
# used within libraries that use the same type tag within multiple
# compile units.
opm_add_test(lens_immiscible_ecfv_ad_mcu
             ONLY_COMPILE
             SOURCES
                tests/lens_immiscible_ecfv_ad_cu1.cc
                tests/lens_immiscible_ecfv_ad_cu2.cc
                tests/lens_immiscible_ecfv_ad_main.cc)

opm_add_test(finger_immiscible_ecfv
             CONDITION ${DUNE_ALUGRID_FOUND})

opm_add_test(finger_immiscible_vcfv
             CONDITION ${DUNE_ALUGRID_FOUND})

opm_add_test(finger_immiscible_ecfv_adaptive
             EXE_NAME finger_immiscible_ecfv
             CONDITION ${DUNE_ALUGRID_FOUND} AND ${DUNE_FEM_FOUND}
             NO_COMPILE
             TEST_ARGS --enable-grid-adaptation=true --end-time=25e3)

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
             powerinjection_forchheimer_fd
             powerinjection_forchheimer_ad
             powerinjection_darcy_fd
             powerinjection_darcy_ad
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
             CONDITION ${DUNE_ALUGRID_FOUND}
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
             CONDITION ${MPI_FOUND}
             DRIVER_ARGS --parallel-simulation=4)

# test for the parallelization of the vertex centered finite volume
# discretization (using BiCGSTAB + ILU0)
opm_add_test(obstacle_immiscible_parallel
             EXE_NAME obstacle_immiscible
             NO_COMPILE
             PROCESSORS 4
             CONDITION ${MPI_FOUND}
             DRIVER_ARGS --parallel-simulation=4
             TEST_ARGS --end-time=1 --initial-time-step-size=1)

# test for the parallel AMG linear solver using the vertex centered
# finite volume discretization
opm_add_test(lens_immiscible_vcfv_fd_parallel
             EXE_NAME lens_immiscible_vcfv_fd
             NO_COMPILE
             PROCESSORS 4
             CONDITION ${MPI_FOUND}
             DRIVER_ARGS --parallel-simulation=4
             TEST_ARGS --end-time=250 --initial-time-step-size=250)

opm_add_test(lens_immiscible_vcfv_ad_parallel
             EXE_NAME lens_immiscible_vcfv_ad
             NO_COMPILE
             PROCESSORS 4
             CONDITION ${MPI_FOUND}
             DRIVER_ARGS --parallel-simulation=4
             TEST_ARGS --end-time=250 --initial-time-step-size=250)

opm_add_test(lens_immiscible_ecfv_ad_parallel
             EXE_NAME lens_immiscible_ecfv_ad
             NO_COMPILE
             PROCESSORS 4
             CONDITION ${MPI_FOUND}
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
             TEST_ARGS --pvs-verbosity=2 --end-time=30000)

opm_add_test(tutorial1
             SOURCES tutorial/tutorial1.cc)

opm_add_test(test_tasklets
             DRIVER_ARGS --plain)
