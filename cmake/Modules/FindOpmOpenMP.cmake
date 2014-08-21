# - Use OpenMP features
#
# Synopsis:
#
#	findOpmOpenMP()
#
# Note: Compiler flags are always added globally, to avoid ABI
# incompatibility problems.
#
# It is assumed that the following variables are available
#
#	OPENMP_QUIET      Verbosity level of the parent's find module
#	OPENMP_LIBRARIES  List of libraries to which OpenMP will be added
#
include (AddOptions)
include (UseCompVer)
is_compiler_gcc_compatible ()

# user code can be conservative by setting USE_OPENMP_DEFAULT
if (NOT DEFINED USE_OPENMP_DEFAULT)
  set(USE_OPENMP_DEFAULT ON)
endif()
option (USE_OPENMP "Enable OpenMP for parallelization" ${USE_OPENMP_DEFAULT})

if (USE_OPENMP)
  # enabling OpenMP is supposedly enough to make the compiler link with
  # the appropriate libraries
  find_package (OpenMP ${OPENMP_QUIET})
  list (APPEND OPMOPENMP_LIBRARIES ${OpenMP_LIBRARIES})
  if (OPENMP_FOUND)
    add_options(C ALL_BUILDS "${OpenMP_C_FLAGS}")
    add_options(CXX ALL_BUILDS "${OpenMP_CXX_FLAGS}")
    set (HAVE_OPENMP 1)
  endif()
else()
  message (STATUS "OpenMP: disabled")
endif()
