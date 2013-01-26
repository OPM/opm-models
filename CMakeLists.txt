##############
# general stuff. we extract the all information possible from the
# dune.module file.
cmake_minimum_required(VERSION 2.8)

# read the dune.module file
file(READ "${CMAKE_SOURCE_DIR}/dune.module" DUNE_MODULE)

# find the project name
string(REGEX REPLACE ".*Module:[ ]*([^\n]+).*" "\\1" ProjectName "${DUNE_MODULE}")

# set the project version. also split the version string into MAJOR.MINOR.REVISON
string(REGEX REPLACE ".*Version:[ ]*([^ \n]+).*" "\\1" ProjectVersion "${DUNE_MODULE}")
string(REGEX REPLACE "^([0-9]*)\\..*\$" "\\1" EWOMS_VERSION_MAJOR "${ProjectVersion}")
string(REGEX REPLACE "^[0-9]*\\.([0-9]*).*\$" "\\1" EWOMS_VERSION_MINOR "${ProjectVersion}")
string(REGEX REPLACE "^[0-9]*\\.[0-9]*\\.([0-9]*).*\$" "\\1" EWOMS_VERSION_REVISION "${ProjectVersion}")

# if the regular expression for the revision did not match, we use "0"
# as the revision number. (we silently assume, that the regexps for
# the major and minor version match.)
if ("${EWOMS_VERSION_REVISION}" STREQUAL "${ProjectVersion}")
  set(EWOMS_VERSION_REVISION "0")
endif("${EWOMS_VERSION_REVISION}" STREQUAL "${ProjectVersion}")

# set the maintainer email (the default Dune autotools build system
# assumes that dune.module's 'Maintainer' field only contains the
# email address of the maintainer. Using the format 'Maintainer:
# Maintainer Name <maintainer@address.org>' makes the DUNE autotools
# build system choke, so we introduce a new field 'MaintainerName'
# which is ignored by the DUNE autotools build system.)
string(REGEX REPLACE ".*MaintainerName:[ ]*([^\n]+).*" "\\1" ProjectMaintainer "${DUNE_MODULE}")
string(REGEX REPLACE ".*Maintainer:[ ]*([^\n]+).*" "\\1" ProjectMaintainerEmail "${DUNE_MODULE}")

# find codename string
string(REGEX REPLACE ".*Codename:[ ]*([^\n]+).*" "\\1" ProjectCodename "${DUNE_MODULE}")

##############

project(${ProjectName} CXX)
# needed for tests like pthread and BLAS
enable_language(C)

##############
# make sure our own modules will be found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# needed for Debian installation scheme
include(UseMultiArch)

# eWoms unit-testing infrastructure
include(EwomsAddTest)
 
##############
# Set CMake policies

# resolve library paths to the policy introduced by CMake 2.6. For
# details, see
# http://www.cmake.org/cmake/help/cmake-2.6.html#policy:CMP0003
cmake_policy(SET CMP0003 NEW)

##############
# Find the required packages

# set required compiler flags for C++11 (former C++0x)
find_package(CXX11Features)
find_package(SharedPtr)

find_package(DUNE_common REQUIRED)
find_package(DUNE_geometry REQUIRED)
find_package(DUNE_grid REQUIRED)
find_package(DUNE_localfunctions REQUIRED)
find_package(DUNE_istl REQUIRED)

##############
# Find the optional packages
find_package(MPI)
find_package(Alberta)
find_package(UG)
find_package(ALUGrid)
find_package(METIS)
find_package(SuperLU)
find_package(Valgrind)

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
    set(${ConfigHName} "${${CMakeName}}")
  else(${CMakeName})
    set(${ConfigHName} 0)
  endif(${CMakeName})
endmacro(EwomsSetConfigHVar)

macro(EwomsSetBoolConfigHVar ConfigHName CMakeName)
  if(${CMakeName})
    set(${ConfigHName} 1)
  else(${CMakeName})
    set(${ConfigHName} 0)
  endif(${CMakeName})
endmacro(EwomsSetBoolConfigHVar)

EwomsSetBoolConfigHVar(HAVE_DUNE_COMMON         DUNE_common_FOUND)
EwomsSetBoolConfigHVar(HAVE_DUNE_GEOMETRY       DUNE_geometry_FOUND)
EwomsSetBoolConfigHVar(HAVE_DUNE_GRID           DUNE_grid_FOUND)
EwomsSetBoolConfigHVar(HAVE_DUNE_LOCALFUNCTIONS DUNE_localfunctions_FOUND)
EwomsSetBoolConfigHVar(HAVE_DUNE_ISTL           DUNE_istl_FOUND)
EwomsSetBoolConfigHVar(HAVE_SUPERLU             SUPERLU_FOUND)
EwomsSetBoolConfigHVar(HAVE_VALGRIND            VALGRIND_INCLUDE_DIR)

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
  set (DUNE_DIR ${CMAKE_INSTALL_LIBDIR_NOARCH}/dunecontrol CACHE PATH "Parent directory of DUNE modules")
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

##############
# deal with MPI
if(MPI_FOUND)
   set(TMP ${MPI_COMPILE_FLAGS})
   separate_arguments(TMP)
   add_definitions(${TMP})
# add_definitions(-DModelP) # tell UG that the model is parallelized
   set(EwomsLinkLibraries ${EwomsLinkLibraries} ${MPI_LIBRARIES})
   set(EwomsIncludeDirectories ${EwomsIncludeDirectories} ${MPI_INCLUDE_PATH})
endif(MPI_FOUND)
EwomsSetBoolConfigHVar(HAVE_MPI MPI_FOUND)
##############

##############
# documentation
configure_file (
        ${PROJECT_SOURCE_DIR}/doc/doxygen/Doxyfile.cmake
        ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile
        @ONLY
        )
find_package (Doxygen)
if (DOXYGEN_FOUND)
        add_custom_target (doc
                COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile
                SOURCES ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc/doxygen
                COMMENT "Generating API documentation with Doxygen"
                VERBATIM
                )
        set (_formats html)
        foreach (format IN LISTS _formats)
                string (TOUPPER ${format} FORMAT)
                install (
                        DIRECTORY ${PROJECT_BINARY_DIR}/doc/doxygen/${format}
                        DESTINATION share/doc/ewoms/
                        COMPONENT ${format}
                        OPTIONAL
                        )
                # target to install just HTML documentation
                add_custom_target (install-${format}
                        COMMAND ${CMAKE_COMMAND} -DCOMPONENT=${format} -P cmake_install.cmake
                        COMMENT Installing ${FORMAT} documentation
                        VERBATIM
                        )
                # since the documentation is optional, it is not automatically built
                add_dependencies (install-${format} doc)
        endforeach (format)
endif (DOXYGEN_FOUND)
##############

##############
# header file installation

# find all the headers. The CMake documentation recommends against
# using globs to enumerate source code, but if you list the files
# explicitly you'll need to change the build files every time you add
# to the project as well as having to rebuild completely...
file (GLOB_RECURSE ewoms_HEADERS "${CMAKE_SOURCE_DIR}/ewoms/[^.]*.hh")

foreach (_hdr IN LISTS ewoms_HEADERS)
	get_filename_component (_dir ${_hdr} PATH)
	file (RELATIVE_PATH _rel_dir "${PROJECT_SOURCE_DIR}" "${_dir}")
	install (
		FILES ${_hdr}
		DESTINATION include/${_rel_dir}
		)
endforeach (_hdr)
install (
	FILES ${CMAKE_SOURCE_DIR}/dune.module
	DESTINATION ${CMAKE_INSTALL_LIBDIR_NOARCH}/dunecontrol/ewoms
	)
message (STATUS "This build defaults to installing in ${CMAKE_INSTALL_PREFIX}")

# use the same name for the source tree as AutoTools so that we can
# share the config file template between the two systems
get_filename_component(abs_top_srcdir "${CMAKE_CURRENT_LIST_FILE}" PATH)

# extract the definitions that various modules have added
get_directory_property (EWOMS_COMP_DEFS COMPILE_DEFINITIONS)
foreach(d ${EWOMS_COMP_DEFS})
  set(ALL_PKG_CPPFLAGS "${ALL_PKG_CPPFLAGS}-D${d} ")
endforeach()

##############
# write the config.h file to disk 
configure_file(config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h )

# make variable substitution and write config mode file
configure_file(EwomsConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/EwomsConfig.cmake @ONLY)

#add_definitions(-DHAVE_CONFIG_H)
##############

##############
# set up CTest 
enable_testing()
include(CTest)

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

##############
# tell cmake that we've got a few subdirectories. (that's the
# directories where the unit tests and tutorials are.)
add_subdirectory("test")
add_subdirectory("tutorial")
##############
