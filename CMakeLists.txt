# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

cmake_minimum_required(VERSION 3.0)

# specify the option that allows to select the build system. if
# dunecontrol is used, we default to the dune build system. if not, we
# use the OPM home-brewn one
set(USE_DUNE_BUILDSYSTEM_DEFAULT OFF)
if (DEFINED ENV{DUNE_CONTROL_PATH})
  set(USE_DUNE_BUILDSYSTEM_DEFAULT ON)
endif()
option(USE_DUNE_BUILDSYSTEM "Use the Dune build system?" ${USE_DUNE_BUILDSYSTEM_DEFAULT})

if (USE_DUNE_BUILDSYSTEM)

######
# use DUNE's build system
######

# set up project and specify the minimum cmake version
project("ewoms" C CXX)

# Sibling build
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
if(SIBLING_SEARCH)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  foreach(MOD dune-common dune-geometry dune-grid dune-localfunctions dune-istl opm-common libecl opm-parser opm-grid opm-core opm-material)
    # do not overwrite the location of a module if it was explicitly
    # specified by the user
    if(${MOD}_DIR)
      continue()
    endif()

    # Try various possible locations for the build directory of the dependency
    foreach(BUILD_DIR "${_leaf_dir_name}" "build-cmake" "build" "../${MOD}-build" "../build-${MOD}" "../build" ".")
      if(EXISTS "${_modules_dir}/${MOD}/${BUILD_DIR}/${MOD}-config.cmake")
        set(${MOD}_DIR ${_modules_dir}/${MOD}/${BUILD_DIR})
        break()
      endif()
    endforeach()
  endforeach()
endif()

if(dune-common_DIR AND NOT IS_DIRECTORY ${dune-common_DIR})
  message(WARNING "Value ${dune-common_DIR} passed to variable"
    " dune-common_DIR is not a directory")
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

# Set CMP0053 (how to handle escape sequences in strings) to the new
# behavior to avoid a pretty annoying cmake warning if a library is
# defined in the toplevel CMakeLists.txt. This should probably be
# considered to be a bug in the dune build system. Note that the old
# behaviour will most likely also work fine, but the result of setting
# this policy to NEW is most likely what's intended.
if (POLICY CMP0053)
  cmake_policy(SET CMP0053 NEW)
endif()

# find the build system (i.e., dune-common) and set cmake's module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH})

# include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# include the additional generic OPM cmake macros
include(OpmMacros)

# recursively mark all header files beneath the "ewoms" and "ebos"
# directories for installation.
opm_recusive_export_all_headers("ewoms" "ebos")

# the cmake modules get a special treatment
opm_export_cmake_modules()

# we want all features detected by the build system to be enabled,
# thank you!
dune_enable_all_packages()

# copy the data files (i.e., grid files, etc) needed to run the unit
# tests from the source to the build directory
opm_recusive_copy_testdata("tests/data/*" "tests/referencesolutions/*")
file(COPY "${CMAKE_SOURCE_DIR}/tests/data/" DESTINATION "${CMAKE_BINARY_DIR}/data")

# extract the code name of the current ewoms revision from the
# corresponding dune.module file
set(ewoms_DUNE_MODULE_FILE "${CMAKE_SOURCE_DIR}/dune.module")
if(NOT EXISTS "${ewoms_DUNE_MODULE_FILE}")
  message(FATAL_ERROR "File '${ewoms_DUNE_MODULE_FILE}' does not exist!")
endif()

extract_line("Codename:" EWOMS_VERSION_CODENAME "${ewoms_DUNE_MODULE_FILE}")
if(EWOMS_VERSION_CODENAME)
  set(EWOMS_VERSION_STRING "${EWOMS_VERSION} (\"${EWOMS_VERSION_CODENAME}\")")
  string(REPLACE "\"" "\\\"" EWOMS_VERSION_STRING "${EWOMS_VERSION_STRING}")
else()
  set(EWOMS_VERSION_STRING "${EWOMS_VERSION}")
endif()

# eWoms comes with its custom test driver script
opm_set_test_driver("${PROJECT_SOURCE_DIR}/bin/runtest.sh" "--simulation")

# the general-purpose ECL simulator ("ebos" == "&ecl &black-&oil &simulator)
opm_add_application(ebos
                    SOURCES ebos/ebos.cc
                    EXE_NAME ebos
                    CONDITION opm-grid_FOUND AND opm-parser_FOUND AND opm-core_FOUND)

if(opm-grid_FOUND AND opm-core_FOUND AND opm-parser_FOUND)
  install(TARGETS ebos DESTINATION bin)
endif()

# the ART to DGF file format conversion utility
opm_add_application(art2dgf
                    SOURCES art2dgf/art2dgf.cc
                    EXE_NAME art2dgf)

add_dependencies(test-suite art2dgf)
install(TARGETS art2dgf DESTINATION bin)

opm_add_test(art2dgf
  NO_COMPILE
  DRIVER_ARGS --plain
  TEST_ARGS "data/fracture-raw.art")

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
             CONDITION dune-alugrid_FOUND)

opm_add_test(finger_immiscible_vcfv
             CONDITION dune-alugrid_FOUND)

opm_add_test(finger_immiscible_ecfv_adaptive
             EXE_NAME finger_immiscible_ecfv
             CONDITION dune-alugrid_FOUND AND dune-fem_FOUND
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

# a test which ensures that black-oil specific fluid state used by the
# black-oil model conforms to the generic fluid state API.
opm_add_test(test_blackoilfluidstate
  DRIVER_ARGS --plain)

opm_add_test(fracture_discretefracture
             CONDITION dune-alugrid_FOUND
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
opm_add_test(lens_immiscible_vcfv_fd_parallel
             EXE_NAME lens_immiscible_vcfv_fd
             NO_COMPILE
             PROCESSORS 4
             CONDITION MPI_FOUND
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
             TEST_ARGS --pvs-verbosity=2 --end-time=30000)

opm_add_test(tutorial1
             SOURCES tutorial/tutorial1.cc)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

else()

######
# use the traditional OPM build system
######

# Mandatory call to project
project(ewoms CXX)

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
      set(opm-common_DIR "${_parent_full_dir}/opm-common}")
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

# the general-purpose ECL simulator ("ebos" == "&ecl &black-&oil &simulator)
EwomsAddApplication(ebos
                    SOURCES ebos/ebos.cc
                    EXE_NAME ebos
                    CONDITION ${OPM_GRID_FOUND} AND ${OPM_PARSER_FOUND} AND ${OPM_CORE_FOUND})

if(OPM_GRID_FOUND AND OPM_CORE_FOUND AND OPM_PARSER_FOUND)
  install(TARGETS ebos DESTINATION bin)
endif()

# the ART to DGF file format conversion utility
EwomsAddApplication(art2dgf
                    SOURCES art2dgf/art2dgf.cc
                    EXE_NAME art2dgf)

add_dependencies(test-suite art2dgf)

opm_add_test(art2dgf
  NO_COMPILE
  DRIVER_ARGS --plain
  TEST_ARGS "data/fracture-raw.art")

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

endif()
