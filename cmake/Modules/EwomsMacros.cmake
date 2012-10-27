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

  # Base path to look for libraries and includes
  if(${EwomsModule}_DIR)
    list(APPEND EwomsModulePath ${${EwomsModule}_DIR})
  endif(${EwomsModule}_DIR)
  if(${EwomsFramework}_DIR)
    list(APPEND EwomsModulePath "${${EwomsFramework}_DIR}/${EwomsModuleName}")
  endif(${EwomsFramework}_DIR)

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
  message("INCLUDE_PATHS: ${EwomsIncludePath}")
  message("${Inc}: ${${Inc}}")
  
  if (NOT ${Inc})
    find_path(${Inc}
              ${HeaderName}
              PATHS ${EwomsIncludePath})
  endif()
  message("${Inc}: ${${Inc}}")

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
  set(${EwomsModule}_FOUND ${EwomsFound} CACHE BOOL INTERNAL)
  set(${EwomsModule}_LIBRARIES ${EwomsLibraries} CACHE STRING INTERNAL)

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
      message(FATAL_ERROR "${EwomsPathMessage}")
    else(EwomsLibsFound)
      message(FATAL_ERROR "${EwomsPathMessage} ${EwomsFailedLibsMessage}")
    endif(EwomsLibsFound)
  endif(NOT EwomsFound AND ${EwomsModule}_FIND_REQUIRED)
endmacro(EwomsCheckFound)
