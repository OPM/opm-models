##############
# general stuff
cmake_minimum_required(VERSION 2.8)
set(ProjectName            "eWoms")
set(ProjectVersion         "2.2-git")
set(ProjectMaintainer      "Andreas Lauser")
set(ProjectMaintainerEmail "Andreas.Lauser_at_iws dot uni-stuttgart dot de")
project(${ProjectName} CXX)
# needed for tests like pthread and BLAS
enable_language(C)

##############
# make sure our own modules will be found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

##############
# Set the policy how CMake resolves library paths to the
# policy introduced by CMake 2.6. For details, see
# http://www.cmake.org/cmake/help/cmake-2.6.html#policy:CMP0003
cmake_policy(SET CMP0003 NEW)

# set required compiler flags for C++11 (former C++0x)
find_package(CXX11Features)
find_package(SharedPtr)

##############
# Find the required packages
find_package(DUNE_grid REQUIRED)
find_package(DUNE_istl REQUIRED)
find_package(DUNE_localfunctions REQUIRED)
find_package(DUNE_geometry REQUIRED)
find_package(DUNE_common REQUIRED)

##############
# Find the optional packages
find_package(MPI)
#find_package(Boost)
find_package(Alberta)
find_package(UG)
find_package(ALUGrid)
find_package(METIS)
find_package(SuperLU)

##############
# use this macros in the CMakelists of the subdirectories. 

# -> for TARGET_LINK_LIBRARIES
set(EwomsLinkLibraries 
    ${DUNE_common_LIBRARIES}
    ${DUNE_geometry_LIBRARIES}
    ${DUNE_grid_LIBRARIES})


# -> for INCLUDE_DIRECTORIES
set(EwomsIncludeDirectories 
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}
    ${DUNE_grid_INCLUDE_DIRS}
    ${DUNE_geometry_INCLUDE_DIRS}
    ${DUNE_common_INCLUDE_DIRS}
    ${DUNE_istl_INCLUDE_DIRS}
    ${DUNE_localfunctions_INCLUDE_DIRS})

if(BOOST_FOUND)
  set(EwomsLinkDirectories "${EwomsLinkDirectories} ${Boost_LIBRARY_DIRS}")
  set(EwomsLinkLibraries "${EwomsLinkLibraries} ${Boost_LIBRARIES}")
  set(EwomsIncludeDirectories "${EwomsIncludeDirectories} ${Boost_INCLUDE_DIR}")
endif(BOOST_FOUND)

if(SUPERLU_FOUND)
  set(EwomsLinkLibraries ${EwomsLinkLibraries} ${SUPERLU_LIBS})
  set(EwomsIncludeDirectories ${EwomsIncludeDirectories} ${SUPERLU_INCLUDE_DIRS})
endif(SUPERLU_FOUND)

##############
# set appropriate compiler flags for debug/release compilation modes
add_definitions("-std=c++0x -Wall -Wno-sign-compare -fno-strict-aliasing")

if(("${CMAKE_BUILD_TYPE}" STREQUAL "release") OR (NOT (DEFINED "CMAKE_BUILD_TYPE")))
  # release mode
  add_definitions("-O3 -march=native")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "ctest")
  # CTest debug mode
  add_definitions(-DDEBUG=1)
  add_definitions("-g -fprofile-arcs -ftest-coverage")
  set(EwomsLinkLibraries "${EwomsLinkLibraries}" "-lm -fprofile-arcs -ftest-coverage")
else()
  # CTest debug mode
  add_definitions(-DDEBUG=1)
  add_definitions("-g")
endif( )

add_definitions(${CUSTOM_CXX_FLAGS}) 

# with the config.h include file...
macro(EwomsSetConfigHVar ConfigHName CMakeName)
  if(${CMakeName})
    set(${ConfigHName} ${${CMakeName}})
  else(${CMakeName})
    set(${ConfigHName} 0)
  endif(${CMakeName})
endmacro(EwomsSetConfigHVar)

macro(EwomsSetBoolConfigHVar ConfigHName CMakeName)
  if(${CMakeName})
    if (${${CMakeName}})
        set(${ConfigHName} 1)
    else (${${CMakeName}})
        set(${ConfigHName} 0)
    endif(${${CMakeName}})
  else(${CMakeName})
    set(${ConfigHName} 0)
  endif(${CMakeName})
endmacro(EwomsSetBoolConfigHVar)

EwomsSetBoolConfigHVar(HAVE_BOOST               Boost_FOUND)
EwomsSetBoolConfigHVar(HAVE_DUNE                DUNE_common_FOUND)
EwomsSetBoolConfigHVar(HAVE_DUNE_GRID           DUNE_grid_FOUND)
EwomsSetBoolConfigHVar(HAVE_DUNE_ISTL           DUNE_istl_FOUND)
EwomsSetBoolConfigHVar(HAVE_DUNE_LOCALFUNCTIONS DUNE_localfunctions_FOUND)
EwomsSetBoolConfigHVar(HAVE_DUNE_PDELAB         DUNE_pdelab_FOUND)
EwomsSetBoolConfigHVar(HAVE_SUPERLU             SUPERLU_FOUND)

EwomsSetConfigHVar(PROJECT_NAME             ProjectName)
EwomsSetConfigHVar(PROJECT_VERSION          ProjectVersion)
EwomsSetConfigHVar(PROJECT_MAINTAINER       ProjectMaintainer)
EwomsSetConfigHVar(PROJECT_MAINTAINER_EMAIL ProjectMaintainerEmail)

# c++ 2011 features
EwomsSetBoolConfigHVar(HAVE_NULLPTR HAVE_NULLPTR)
EwomsSetBoolConfigHVar(HAVE_ARRAY HAVE_ARRAY)
EwomsSetBoolConfigHVar(HAVE_ATTRIBUTE_ALWAYS_INLINE HAVE_ATTRIBUTE_ALWAYS_INLINE)
EwomsSetBoolConfigHVar(HAS_ATTRIBUTE_UNUSED HAS_ATTRIBUTE_UNUSED)
EwomsSetBoolConfigHVar(HAS_ATTRIBUTE_DEPRECATED HAS_ATTRIBUTE_DEPRECATED)
EwomsSetBoolConfigHVar(HAS_ATTRIBUTE_DEPRECATED_MSG HAS_ATTRIBUTE_DEPRECATED_MSG)
EwomsSetBoolConfigHVar(HAVE_INTEGRAL_CONSTANT HAVE_INTEGRAL_CONSTANT)
EwomsSetBoolConfigHVar(HAVE_STATIC_ASSERT HAVE_STATIC_ASSERT)
EwomsSetBoolConfigHVar(HAVE_VARIADIC_TEMPLATES HAVE_VARIADIC_TEMPLATES)
EwomsSetBoolConfigHVar(HAVE_VARIADIC_CONSTRUCTOR_SFINAE HAVE_VARIADIC_CONSTRUCTOR_SFINAE)
EwomsSetBoolConfigHVar(HAVE_RVALUE_REFERENCES HAVE_RVALUE_REFERENCES)
EwomsSetBoolConfigHVar(HAVE_TUPLE HAVE_TUPLE)
EwomsSetBoolConfigHVar(HAVE_TR1_TUPLE HAVE_TR1_TUPLE)

##############
# default directory to find control files if installed by package
if (UNIX)
  set (DUNE_DIR "/usr/lib/dunecontrol" CACHE PATH "Parent directory of DUNE modules")
endif (UNIX)

# add dune-common version from dune.module to config.h
if(DUNE_common_DIR)
  file(READ "${DUNE_common_DIR}/dune.module" DUNE_COMMON_MODULE)
else()
  # note that forward slash works on WIN32 too
  if(NOT EXISTS "${DUNE_DIR}/dune-common/dune.module")
    message(FATAL_ERROR "No dune-common directory found")
  endif(NOT EXISTS "${DUNE_DIR}/dune-common/dune.module")
  file(READ "${DUNE_DIR}/dune-common/dune.module" DUNE_COMMON_MODULE)
endif(DUNE_common_DIR)

# find version string
string(REGEX REPLACE ".*Version:[ ]*([^ \n]+).*" "\\1" DUNE_COMMON_VERSION "${DUNE_COMMON_MODULE}")
string(REGEX REPLACE "([0-9]).*" "\\1" DUNE_COMMON_VERSION_MAJOR "${DUNE_COMMON_VERSION}")
string(REGEX REPLACE "[0-9]*\\.([0-9]).*" "\\1" DUNE_COMMON_VERSION_MINOR "${DUNE_COMMON_VERSION}")
string(REGEX REPLACE "[0-9]*\\.[0-9]*\\.([0-9]).*" "\\1" DUNE_COMMON_VERSION_REVISION "${DUNE_COMMON_VERSION}")

# remove false matches
string(REGEX MATCH "[^0-9]" NON_NUMBER_CHARACTER "${DUNE_COMMON_VERSION_MINOR}")
if(NON_NUMBER_CHARACTER)
  set(DUNE_COMMON_VERSION_MINOR "0")
endif(NON_NUMBER_CHARACTER)
string(REGEX MATCH "[^0-9]" NON_NUMBER_CHARACTER "${DUNE_COMMON_VERSION_REVISION}")
if(NON_NUMBER_CHARACTER)
  set(DUNE_COMMON_VERSION_REVISION "0")
endif(NON_NUMBER_CHARACTER)
##############

##############
# add eWoms version from dune.module to config.h
file(READ "${CMAKE_SOURCE_DIR}/dune.module" EWOMS_MODULE)

# find version string
string(REGEX REPLACE ".*Version:[ ]*([^ \n]+).*" "\\1" EWOMS_VERSION "${EWOMS_MODULE}")
string(REGEX REPLACE "([0-9]).*" "\\1" EWOMS_VERSION_MAJOR "${EWOMS_VERSION}")
string(REGEX REPLACE "[0-9]*\\.([0-9]).*" "\\1" EWOMS_VERSION_MINOR "${EWOMS_VERSION}")
string(REGEX REPLACE "[0-9]*\\.[0-9]*\\.([0-9]).*" "\\1" EWOMS_VERSION_REVISION "${EWOMS_VERSION}")

# remove false matches
string(REGEX MATCH "[^0-9]" NON_NUMBER_CHARACTER "${EWOMS_VERSION_MINOR}")
if(NON_NUMBER_CHARACTER)
  set(EWOMS_VERSION_MINOR "0")
endif(NON_NUMBER_CHARACTER)
string(REGEX MATCH "[^0-9]" NON_NUMBER_CHARACTER "${EWOMS_VERSION_REVISION}")
if(NON_NUMBER_CHARACTER)
  set(EWOMS_VERSION_REVISION "0")
endif(NON_NUMBER_CHARACTER)

# find codename string
string(REGEX REPLACE ".*Codename:[ ]*([^\n]+).*" "\\1" EWOMS_CODENAME "${EWOMS_MODULE}")
set(EWOMS_CODENAME "${EWOMS_CODENAME}")

##############

##############
# adapt build system to detected packages

# deal with UG
if(UG_FOUND)
  set(EwomsIncludeDirectories ${EwomsIncludeDirectories} ${UG_INCLUDE_DIRS})
  set(EwomsLinkLibraries ${EwomsLinkLibraries} ${UG_LIBRARIES})
endif(UG_FOUND)
EwomsSetBoolConfigHVar(HAVE_UG UG_FOUND)

# deal with ALUGrid
if(ALUGRID_FOUND)
  set(EwomsIncludeDirectories ${EwomsIncludeDirectories} ${ALUGRID_INCLUDES})
  set(EwomsLinkLibraries ${EwomsLinkLibraries} ${ALUGRID_LIB})
endif()
EwomsSetBoolConfigHVar(HAVE_ALUGRID ALUGRID_FOUND)

# deal with Alberta
if(Alberta_FOUND)
 set(EwomsIncludeDirectories ${EwomsIncludeDirectories} ${Alberta_INCLUDE_DIRS})
 set(EwomsLinkLibraries ${EwomsLinkLibraries} ${Alberta_LIBRARIES})
endif(Alberta_FOUND)
EwomsSetBoolConfigHVar(HAVE_ALBERTA Alberta_FOUND)

# deal with METIS
if(METIS_FOUND)
 set(EwomsIncludeDirectories ${EwomsIncludeDirectories} ${METIS_INCLUDE_DIRS})
 set(EwomsLinkLibraries ${EwomsLinkLibraries} ${METIS_LIBRARIES})
endif(METIS_FOUND)
EwomsSetBoolConfigHVar(HAVE_METIS METIS_FOUND)
##############

if(MPI_FOUND)
   set(TMP ${MPI_COMPILE_FLAGS})
   separate_arguments(TMP)
   add_definitions(${TMP})
# add_definitions(-DModelP) # tell UG that the model is parallelized
   set(EwomsLinkLibraries ${EwomsLinkLibraries} ${MPI_LIBRARIES})
   set(EwomsIncludeDirectories ${EwomsIncludeDirectories} ${MPI_INCLUDE_PATH})
endif(MPI_FOUND)
EwomsSetBoolConfigHVar(HAVE_MPI MPI_FOUND)

# actually write the config.h file to disk 
configure_file(config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h )
#add_definitions(-DHAVE_CONFIG_H)

##############
# tell cmake that we've got a few subdirectories. (that's the
# directories where the actual programs are)
add_subdirectory("test")
add_subdirectory("tutorial")

# copy the testing script
make_directory(bin)
make_directory(referencesolutions)
make_directory(grids)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/bin/runtest.sh
           DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/bin/fuzzycomparevtu.py
           DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin)
file(GLOB vtus
           "${CMAKE_CURRENT_SOURCE_DIR}/test/referencesolutions/*.vtu" 
           "${CMAKE_CURRENT_SOURCE_DIR/}test/referencesolutions/*.vtp")
foreach(file ${vtus})
    file(COPY ${file} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/referencesolutions)
endforeach(file)

file(GLOB grids "test/implicit/grids/*.dgf" "test/decoupled/*/grids/*.dgf" "test/implicit/grids/*.art")
foreach(file ${grids})
    file(COPY ${file} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/grids)
endforeach(file)

# set up CTest 
enable_testing()
include(CTest)

# we add the water-air test first because it takes longest by far
add_test(waterair_pvs_ni bin/runtest.sh --simulation waterair_pvs_ni)

add_test(test_propertysystem bin/runtest.sh --plain test_propertysystem)
add_test(test_spline bin/runtest.sh --plain test_spline)
add_test(test_fluidsystems bin/runtest.sh --plain test_fluidsystems)
add_test(test_immiscibleflash bin/runtest.sh --plain test_immiscibleflash)
add_test(test_ncpflash bin/runtest.sh --plain test_ncpflash)
add_test(test_pengrobinson bin/runtest.sh --plain test_pengrobinson)
add_test(test_tabulation bin/runtest.sh --plain test_tabulation)

add_test(co2injection_immiscible bin/runtest.sh --simulation co2injection_immiscible)
add_test(co2injection_immiscible_ni bin/runtest.sh --simulation co2injection_immiscible_ni)
add_test(co2injection_flash bin/runtest.sh --simulation co2injection_flash)
add_test(co2injection_flash_ni bin/runtest.sh --simulation co2injection_flash_ni)
add_test(co2injection_ncp bin/runtest.sh --simulation co2injection_ncp)
add_test(co2injection_ncp_ni bin/runtest.sh --simulation co2injection_ncp_ni)
add_test(co2injection_pvs bin/runtest.sh --simulation co2injection_pvs)
add_test(co2injection_pvs_ni bin/runtest.sh --simulation co2injection_pvs_ni)
add_test(cuvette_pvs bin/runtest.sh --simulation cuvette_pvs)
add_test(diffusion_flash bin/runtest.sh --simulation diffusion_flash)
add_test(diffusion_ncp bin/runtest.sh --simulation diffusion_ncp)
add_test(diffusion_pvs bin/runtest.sh --simulation diffusion_pvs)
add_test(finger_immiscible bin/runtest.sh --simulation finger_immiscible)
add_test(fracture_discretefracture bin/runtest.sh --simulation fracture_discretefracture --end-time=1000)
add_test(groundwater_immiscible bin/runtest.sh --simulation groundwater_immiscible --linear-solver-verbosity=2)
add_test(infiltration_pvs bin/runtest.sh --simulation infiltration_pvs)
add_test(lens_immiscible bin/runtest.sh --simulation lens_immiscible)
add_test(lens_richards bin/runtest.sh --simulation lens_richards --newton-write-convergence=true)
add_test(obstacle_immiscible bin/runtest.sh --simulation obstacle_immiscible)
add_test(obstacle_ncp bin/runtest.sh --simulation obstacle_ncp)
add_test(obstacle_pvs bin/runtest.sh --simulation obstacle_pvs)
add_test(outflow_pvs bin/runtest.sh --simulation outflow_pvs)
add_test(powerinjection_darcy bin/runtest.sh --simulation powerinjection_darcy)
add_test(powerinjection_forchheimer bin/runtest.sh --simulation powerinjection_forchheimer)
add_test(reservoir_blackoil bin/runtest.sh --simulation reservoir_blackoil)
add_test(test_quadrature bin/runtest.sh --plain test_quadrature)

add_test(test_diffusion bin/runtest.sh --simulation-diffusion test_diffusion 3)
add_test(test_1p bin/runtest.sh --plain test_1p)
add_test(test_transport bin/runtest.sh --simulation test_transport)
add_test(test_impes bin/runtest.sh --simulation test_impes)
add_test(test_multiphysics2p2c bin/runtest.sh --simulation test_multiphysics2p2c --end-time=2000)

add_test(tutorial_coupled bin/runtest.sh --simulation tutorial_coupled)
add_test(tutorial_decoupled bin/runtest.sh --simulation tutorial_decoupled)

# test the restart capability
add_test(obstacle_pvs_restart bin/runtest.sh --restart obstacle_pvs --pvs-verbosity=2)

# test the parameter system
add_test(obstacle_immiscible_parameters bin/runtest.sh --parameters obstacle_immiscible)

# conditionally add the tests that depend on superLU and ALUGrid
if(SUPERLU_FOUND AND ALUGRID_FOUND)
  add_test(test_navierstokes bin/runtest.sh --simulation test_navierstokes)
endif()

# conditionally add the tests that depend on ALUGrid
if(ALUGRID_FOUND)
  add_test(test_mpfao2p bin/runtest.sh --simulation test_mpfa2p -ModelType MPFAO)
  add_test(test_mpfal2p bin/runtest.sh --simulation test_mpfa2p -ModelType MPFAL)
  add_test(test_mpfal2padaptive  bin/runtest.sh --simulation test_mpfa2p -ModelType MPFALAdaptive --end-time=2e4)
  add_test(test_adaptive2p2c bin/runtest.sh --simulation test_adaptive2p2c)
  add_test(test_impesadaptive bin/runtest.sh --simulation test_impesadaptive)
endif()

# run a test in parallel
if(MPI_FOUND)
  add_test(obstacle_immiscible_parallel bin/runtest.sh --parallel-simulation obstacle_immiscible --end-time=1 --initial-time-step-size=1)
endif()

# conditionally add the tests that depend on superLU
if(SUPERLU_FOUND)
 add_test(test_stokes bin/runtest.sh --simulation test_stokes)
 add_test(test_stokes2c bin/runtest.sh --simulation test_stokes2c)
 add_test(test_stokesni bin/runtest.sh --simulation test_stokesni)
endif()

# use the same name for the source tree as AutoTools so that we can
# share the config file template between the two systems
get_filename_component(abs_top_srcdir "${CMAKE_CURRENT_LIST_FILE}" PATH)

# extract the definitions that various modules have added
get_directory_property (EWOMS_COMP_DEFS COMPILE_DEFINITIONS)
foreach(d ${EWOMS_COMP_DEFS})
  set(ALL_PKG_CPPFLAGS "${ALL_PKG_CPPFLAGS}-D${d} ")
endforeach()

# make variable substitution and write config mode file
configure_file(EwomsConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/EwomsConfig.cmake @ONLY)
