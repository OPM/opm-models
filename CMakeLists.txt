# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# set up project and specify the minimum cmake version
project("ewoms" C CXX)

cmake_minimum_required(VERSION 3.0)

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
