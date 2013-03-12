# helper macro to retrieve a single field of a dune.module file
macro(EwomsGetDuneModuleDirective_ FieldName OutputVariable DuneModuleContents)
  string(REGEX MATCH ".*${FieldName}:[ ]*([^\n]+).*" ${OutputVariable} "${DuneModuleContents}")
  string(REGEX REPLACE ".*${FieldName}:[ ]*([^\n]+).*" "\\1" "${OutputVariable}" "${${OutputVariable}}")
endmacro()

# This macro parses the dune.module file of a Dune module.
#
# Usage:
#
# EwomsParseDuneModuleInfo(DuneModuleName
#              [FILE_NAME "PathTo/dune.module"])
#
# This macro sets the following variables:
#
# ${MODULE_NAME}_NAME: The name of the module which ought to be used when referencing it in texts
# ${MODULE_NAME}_DESCRIPTION: A textual description of the module
# ${MODULE_NAME}_URL: The URL of the module's website
# ${MODULE_NAME}_MAINTAINER_NAME: The real name of the module maintainer(s)
# ${MODULE_NAME}_MAINTAINER_EMAIL: The e-mail address of the module maintainer(s)
# ${MODULE_NAME}_VERSION: version string of the dune module
# ${MODULE_NAME}_VERSION_MAJOR: major version of the dune module
# ${MODULE_NAME}_VERSION_MINOR: minor version of the dune module
# ${MODULE_NAME}_VERSION_REVISION: revision of the dune module
# ${MODULE_NAME}_CODENAME: code name for the version of the module
macro(EwomsParseDuneModuleInfo ModuleName)
  CMAKE_PARSE_ARGUMENTS(
    CURMOD # prefix
    "" # flags
    "FILE_NAME" # one value args
    "" # multi-value args
    ${ARGN})

  # create an uppercase, underscore version of the given module name
  string(TOUPPER "${ModuleName}" MODULE_NAME)
  string(REPLACE "-" "_" MODULE_NAME "${MODULE_NAME}")
  
  if (NOT CURMOD_FILE_NAME)
    set(CURMOD_FILE_NAME "${${MODULE_NAME}_DIR}/dune.module")
  endif()

  # read the dune.module file
  file(READ "${CURMOD_FILE_NAME}" DUNE_MODULE)

  # set the module name
  EwomsGetDuneModuleDirective_("Module" ${MODULE_NAME}_NAME "${DUNE_MODULE}")

  # set the module description
  EwomsGetDuneModuleDirective_("Description" ${MODULE_NAME}__DESCRIPTION "${DUNE_MODULE}")

  # set the URL of the module's website
  EwomsGetDuneModuleDirective_("Url" ${MODULE_NAME}_URL "${DUNE_MODULE}")

  # set the project version. also split the version string into MAJOR.MINOR.REVISON
  EwomsGetDuneModuleDirective_("Version" ${MODULE_NAME}_VERSION "${DUNE_MODULE}")

  string(REGEX REPLACE "^([0-9]*)\\..*\$" "\\1" ${MODULE_NAME}_VERSION_MAJOR "${${MODULE_NAME}_VERSION}")
  string(REGEX REPLACE "^[0-9]*\\.([0-9]*).*\$" "\\1" ${MODULE_NAME}_VERSION_MINOR "${${MODULE_NAME}_VERSION}")
  string(REGEX REPLACE "^[0-9]*\\.[0-9]*\\.([0-9]*).*\$" "\\1" ${MODULE_NAME}_VERSION_REVISION "${${MODULE_NAME}_VERSION}")

  # if the regular expression for the revision did not match, we use "0"
  # as the revision number. (we silently assume, that the regexps for
  # the major and minor version match.)
  if ("${${MODULE_NAME}_VERSION_REVISION}" STREQUAL "${${MODULE_NAME}_VERSION}")
    set(${MODULE_NAME}_VERSION_REVISION "0")
  endif()

  # set the maintainer email (the default Dune autotools build system
  # assumes that dune.module's 'Maintainer' field only contains the
  # email address of the maintainer. Using the format 'Maintainer:
  # Maintainer Name <maintainer@address.org>' makes the DUNE autotools
  # build system choke, so we introduce a new field 'MaintainerName'
  # which is ignored by the DUNE autotools build system.)
  EwomsGetDuneModuleDirective_("MaintainerName" ${MODULE_NAME}_MAINTAINER_NAME "${DUNE_MODULE}")
  EwomsGetDuneModuleDirective_("Maintainer" ${MODULE_NAME}_MAINTAINER_EMAIL "${DUNE_MODULE}")

  # find codename string
  EwomsGetDuneModuleDirective_("Codename" ${MODULE_NAME}_CODENAME "${DUNE_MODULE}")

endmacro()

#############################################################
# This sets up the EwomsMacros for the current CMake module.
# Call EwomsSweep at the end of your CMake module in order to
# not pullute the namespace with unused variables
#############################################################
macro(EwomsSetup 
        CMakeModuleName
        ModuleName
        Framework)
  #############
  # Set some internal variables which are used within 
  # the current CMake Module
  set(EwomsModule        ${CMakeModuleName})
  set(EwomsModuleName    ${ModuleName})
  set(EwomsFramework     ${Framework})

  set(EwomsLibsFound 1)
  set(EwomsLibraryNames)
  set(EwomsFound     0)

  set(EwomsPathMessage 
"Set the ${EwomsModule}_DIR cmake cache entry to the directory 
where the ${EwomsModuleName} libraries reside. Alternatively you can set
the ${EwomsFramework}_DIR entry where all ${EwomsFramework} sub-modules have been compiled.")

  # Base path to look for libraries and includes and auxiliary files
  if(NOT ${EwomsModule}_DIR AND ${EwomsFramework}_DIR)
    set(${EwomsModule}_DIR "${${EwomsFramework}_DIR}/${EwomsModuleName}")
  endif()
  if(${EwomsModule}_DIR)
    list(APPEND EwomsModulePath ${${EwomsModule}_DIR})
  endif()

  # Path to look for includes (->EwomsIncludePath) and libraries (-> EwomsLibraryPath)
  set(EwomsIncludePath "")
  foreach(tmp ${EwomsModulePath})
    list(APPEND EwomsIncludePath "${tmp}" "${tmp}/include")
    list(APPEND EwomsLibraryPath "${tmp}" "${tmp}/lib" "${tmp}/lib64")
  endforeach(tmp)

  set(EwomsLibraries)
  set(EwomsFailedLibraries)
endmacro(EwomsSetup)

#############################################################
# This adds some additional paths to the location where 
# includes and libraries are searched
#############################################################
macro(EwomsAddPathSuffixes 
        IncludeSuffixes
        LibSuffixes)
  foreach(tmp ${EwomsModulePath})
    # deal with the user defined library locations
    foreach(foo ${LibSuffixes})
      list(APPEND EwomsLibraryPath "${tmp}/${foo}")
    endforeach(foo)

    # deal with the user defined include locations
    foreach(foo ${IncludeSuffixes})
      list(APPEND EwomsIncludePath "${tmp}/${foo}")
    endforeach(foo)
  endforeach(tmp)
endmacro(EwomsAddPathSuffixes)
#############################################################
# Find a given library using some reasonable default 
# search paths. Sets Ewoms${LibName}_LIBRARY to the location
# where the library was found and extends the EwomsLibraries
# variable.
#############################################################
macro(EwomsFindLibrary LibName)
  set(Lib ${EwomsModule}_${LibName}_LIBRARY)

  find_library(${Lib}
               ${LibName}
               PATHS ${EwomsLibraryPath}
               PATH_SUFFIXES ".libs"  "lib" "lib32" "lib64"
               NO_DEFAULT_PATH)
  if(NOT ${Lib})
    find_library(${Lib}
                 ${LibName}
                 PATHS ${EwomsLibraryPath}
                 PATH_SUFFIXES ".libs" "lib" "lib32" "lib64")
  endif()

  # unset the ${Lib} variable. since there does not seem to be a way
  # to only remove it from the cache without destroying it in the
  # local scope, we restore it immediately afterwards.
  set(TMP ${${Lib}})
  unset(${Lib} CACHE)
  set(${Lib} ${TMP})

  if(${Lib})
    list(APPEND EwomsLibraries ${${Lib}})
    list(APPEND EwomsLibraryNames ${LibName})
  else(${Lib})
    list(APPEND EwomsFailedLibraries ${LibName})
  endif(${Lib})
endmacro(EwomsFindLibrary)

#############################################################
# Find a given header file using some reasonable default 
# search paths.
#############################################################
macro(EwomsFindIncludeDir HeaderName)
  set(Inc ${EwomsModule}_INCLUDE_DIR)
  find_path(${Inc}
            ${HeaderName}
            PATHS ${EwomsIncludePath} NO_DEFAULT_PATH)
  
  if (NOT ${Inc})
    find_path(${Inc}
              ${HeaderName}
              PATHS ${EwomsIncludePath})
  endif()

  if(${Inc})
    list(APPEND ${EwomsModule}_INCLUDE_DIRS "${${Inc}}")
    list(APPEND EwomsIncludes ${Inc})
  else(${Inc})
    list(APPEND EwomsFailedIncludes ${HeaderName})
  endif(${Inc})
endmacro(EwomsFindIncludeDir)

macro(EwomsFindIncludeBaseDir HeaderName DirSuffix)
  set(Inc ${EwomsModule}_INCLUDE_DIR)
  find_path(${Inc}
            ${HeaderName}
            PATHS ${EwomsIncludePath} NO_DEFAULT_PATH)
  if (NOT ${Inc})
    find_path(${Inc}
              ${HeaderName}
              PATHS ${EwomsIncludePath})
  endif()

  if(${Inc})
    list(APPEND ${EwomsModule}_INCLUDE_DIRS "${${Inc}}/${DirSuffix}")
    list(APPEND EwomsIncludes ${Inc})
  else(${Inc})
    list(APPEND EwomsFailedIncludes ${HeaderName})
  endif(${Inc})
endmacro(EwomsFindIncludeBaseDir)

#############################################################
# Make sure the required libraries were found
#############################################################
macro(EwomsRequiredLibsFound)
  set(EwomsLibsFound 1)
  set(EwomsFailedLibsMessage "Could not find the required libraries ")
  foreach(curLib ${ARGN})
    set(curLibFound 0)
    foreach(tmp ${EwomsLibraryNames})
      if (tmp STREQUAL ${curLib})
        set(curLibFound 1)
      endif (tmp STREQUAL ${curLib})
    endforeach(tmp)
    
    if (NOT curLibFound)
      set(EwomsLibsFound 0)
      set(EwomsFailedLibsMessage "${EwomsFailedLibsMessage} '${curLib}'")
    endif(NOT curLibFound)
  endforeach(curLib)
endmacro(EwomsRequiredLibsFound)

#############################################################
# Make sure the required libraries were found
#############################################################
macro(EwomsIncludeDirsFound)
endmacro(EwomsIncludeDirsFound)

#############################################################
# Make sure everything required was found
#############################################################
macro(EwomsCheckFound)
  # Set the global macros
  set(EwomsFound 0)

  if(EwomsLibsFound AND ${EwomsModule}_INCLUDE_DIR)
    set(EwomsFound 1)
    # log result
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
      "Determing location of ${EwomsModule} succeded:\n"
      "Include directory: ${${EwomsModule}_INCLUDE_DIR}\n"
      "Library directory: ${EwomsLibraries}\n\n")
  else()
    # log errornous result
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
      "Determing location of ${EwomsModule} failed:\n"
      "Include directory: ${${EwomsModule}_INCLUDE_DIR}\n"
      "Library directory: ${EwomsLibraries}\n\n")
  endif(EwomsLibsFound AND ${EwomsModule}_INCLUDE_DIR)
  set(${EwomsModule}_FOUND ${EwomsFound})
  set(${EwomsModule}_LIBRARIES ${EwomsLibraries} CACHE STRING INTERNAL)
  mark_as_advanced(${EwomsModule}_LIBRARIES)
  mark_as_advanced(${EwomsModule}_INCLUDE_DIR)

  # print status message if requested
  if(NOT ${EwomsModule}_FIND_QUIETLY)
    if(EwomsFound)
      message(STATUS "Found ${EwomsModule}")
    else()
      message(STATUS "Could not find ${EwomsModule}")
    endif(EwomsFound)
  endif(NOT ${EwomsModule}_FIND_QUIETLY)

  if(NOT EwomsFound AND ${EwomsModule}_FIND_REQUIRED)
    if(EwomsLibsFound)
      message(FATAL_ERRR "${EwomsPathMessage}")
    else(EwomsLibsFound)
      message(FATAL_ERROR "${EwomsPathMessage} ${EwomsFailedLibsMessage}")
    endif(EwomsLibsFound)
  endif(NOT EwomsFound AND ${EwomsModule}_FIND_REQUIRED)
endmacro(EwomsCheckFound)
