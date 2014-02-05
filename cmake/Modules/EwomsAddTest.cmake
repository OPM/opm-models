# add a unit test. the test is build default, except if the
# DISABLE_TESTING variable is set.
#
# EwomsAddTest(TestName
#              [NO_COMPILE]
#              [ONLY_COMPILE]
#              [EXE_NAME TestExecutableName]
#              [CONDITION ConditionalExpression]
#              [DRIVER_ARGS TestDriverScriptArguments]
#              [SOURCES SourceFile1 SourceFile2 ...])
include(CMakeParseArguments)

macro(EwomsAddTest TestName)
  CMAKE_PARSE_ARGUMENTS(CURTEST
    "NO_COMPILE;ONLY_COMPILE" # flags
    "EXE_NAME" # one value args
    "CONDITION;DRIVER_ARGS;SOURCES" # multi-value args
    ${ARGN})

  set(BUILD_TESTING "${BUILD_TESTING}")

  if (BUILD_TESTING)
    if (NOT CURTEST_EXE_NAME)
      set(CURTEST_EXE_NAME ${TestName})
    endif()
    if (NOT CURTEST_SOURCES)
      set(CURTEST_SOURCES ${CURTEST_EXE_NAME}.cc)
    endif()
    if (NOT CURTEST_DRIVER_ARGS)
      set(CURTEST_DRIVER_ARGS --simulation ${CURTEST_EXE_NAME})
    endif()
    
    if (NOT DISABLE_TESTS)  
      if (NOT CURTEST_ONLY_COMPILE)
        set(SKIP_CUR_TEST "1")
        if (CURTEST_CONDITION STREQUAL "")
          set(SKIP_CUR_TEST "0")
        elseif(${CURTEST_CONDITION})
          set(SKIP_CUR_TEST "0")
        endif()

        if (NOT SKIP_CUR_TEST)
          if (NOT CURTEST_NO_COMPILE)
            add_executable("${CURTEST_EXE_NAME}" ${CURTEST_EXCLUDE_FROM_ALL} ${CURTEST_SOURCES})
            target_link_libraries (${CURTEST_EXE_NAME} ${${project}_LIBRARIES})
          endif()

          add_test(NAME ${TestName} 
            WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
            COMMAND ${CMAKE_SOURCE_DIR}/bin/runtest.sh ${CURTEST_DRIVER_ARGS})
        else ()
          add_test(${TestName} skip_test_dummy)
        endif()
      endif()
    endif()
  endif() # BUILD_TESTING
endmacro()
