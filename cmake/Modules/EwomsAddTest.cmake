# add a unit test. the test is only build by default if the
# ENABLE_TESTING variable was set. syntay
#
# EwomsAddTest(TestName
#              [NO_COMPILE]
#              [EXE_NAME TestExecutableName]
#              [CONDITION ConditionalExpression]
#              [DRIVER_ARGS TestDriverScriptArguments]
#              [SOURCES SourceFile1 SourceFile2 ...])
macro(EwomsAddTest TestName)
  CMAKE_PARSE_ARGUMENTS(CURTEST
    "NO_COMPILE" # flags
    "EXE_NAME" # one value args
    "CONDITION;DRIVER_ARGS;SOURCES" # multi-value args
    ${ARGN})

  if (NOT CURTEST_EXE_NAME)
    set(CURTEST_EXE_NAME ${TestName})
  endif()
  if (NOT CURTEST_SOURCES)
    set(CURTEST_SOURCES ${CURTEST_EXE_NAME}.cc)
  endif()
  if (NOT CURTEST_DRIVER_ARGS)
    set(CURTEST_DRIVER_ARGS --simulation ${CURTEST_EXE_NAME})
  endif()
  
  if (NOT CURTEST_NO_COMPILE)
    if(NOT ENABLE_TESTS)
      # do not include the test binary in 'make all' if
      # EWOMS_ENABLE_TESTS was not specified by the user
      set(CURTEST_EXCLUDE_FROM_ALL "EXCLUDE_FROM_ALL")
    endif()
    add_executable("${CURTEST_EXE_NAME}" ${CURTEST_EXCLUDE_FROM_ALL} ${CURTEST_SOURCES})

    # add required libraries and includes to the build flags 
    target_link_libraries("${CURTEST_EXE_NAME}" ${EwomsLinkLibraries})
    link_directories(${EwomsLinkDirectories})
    include_directories(${EwomsIncludeDirectories})
  endif()
  
  if (ENABLE_TESTS)
    set(SKIP_CUR_TEST "1")
    if ("${CURTEST_CONDITION}" STREQUAL "")
      set(SKIP_CUR_TEST "0")
    else()
      if (${CURTEST_CONDITION})
        set(SKIP_CUR_TEST "0")
      endif()
    endif()
    if (NOT SKIP_CUR_TEST)
      add_test(NAME ${TestName} 
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
        COMMAND bin/runtest.sh ${CURTEST_DRIVER_ARGS})
    else ()
      add_test(${TestName} skip_test_dummy)
    endif()
  endif()
endmacro()
