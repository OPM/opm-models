# add a unit test. the test is build default, except if the
# DISABLE_TESTING variable is set.
#
# EwomsAddTest(TestName
#              [NO_COMPILE]
#              [ONLY_COMPILE]
#              [EXE_NAME TestExecutableName]
#              [CONDITION ConditionalExpression]
#              [DRIVER_ARGS TestDriverScriptArguments]
#              [SOURCES SourceFile1 SourceFile2 ...]
#              [PROCESSORS NumberOfRequiredCores]
#              [DEPENDS TestName1 TestName2 ...])
include(CMakeParseArguments)

macro(EwomsAddTest TestName)
  CMAKE_PARSE_ARGUMENTS(CURTEST
    "NO_COMPILE;ONLY_COMPILE" # flags
    "EXE_NAME;PROCESSORS" # one value args
    "CONDITION;DEPENDS;DRIVER_ARGS;SOURCES" # multi-value args
    ${ARGN})

  set(BUILD_TESTING "${BUILD_TESTING}")

  if (NOT CURTEST_EXE_NAME)
    set(CURTEST_EXE_NAME ${TestName})
  endif()
  if (NOT CURTEST_SOURCES)
    set(CURTEST_SOURCES ${CURTEST_EXE_NAME}.cc)
  endif()
  if (NOT CURTEST_DRIVER_ARGS)
    set(CURTEST_DRIVER_ARGS --simulation ${CURTEST_EXE_NAME})
  endif()
  if (NOT BUILD_TESTING)
    # don't build the tests by _default_ (i.e., when typing
    # 'make'). Though they can still be build using 'make ctests' and
    # they can be build and run using 'make check'
    set(CURTEST_EXCLUDE_FROM_ALL "EXCLUDE_FROM_ALL")
  endif()

  set(SKIP_CUR_TEST "1")
  # the "AND OR " is a hack which is required to prevent CMake from
  # evaluating the condition in the string. (which might
  # evaluate to an empty string even though "${CURTEST_CONDITION}"
  # is not empty.)
  if ("AND OR ${CURTEST_CONDITION}" STREQUAL "AND OR ")
    set(SKIP_CUR_TEST "0")
  elseif(${CURTEST_CONDITION})
    set(SKIP_CUR_TEST "0")
  endif()

  if (NOT SKIP_CUR_TEST)
    if (CURTEST_ONLY_COMPILE)
      add_executable("${CURTEST_EXE_NAME}" ${CURTEST_EXCLUDE_FROM_ALL} ${CURTEST_SOURCES})
      target_link_libraries (${CURTEST_EXE_NAME} ${${project}_LIBRARIES})
    else()
      if (NOT CURTEST_NO_COMPILE)
        add_executable("${CURTEST_EXE_NAME}" ${CURTEST_EXCLUDE_FROM_ALL} ${CURTEST_SOURCES})
        target_link_libraries (${CURTEST_EXE_NAME} ${${project}_LIBRARIES})
        add_dependencies("check" "${CURTEST_EXE_NAME}")

        if(NOT TARGET "ctests")
          add_custom_target("ctests")
        endif()
        add_dependencies("ctests" "${CURTEST_EXE_NAME}")
      endif()

      add_test(NAME ${TestName}
        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
        ${DEPENDS_ON}
        COMMAND ${CMAKE_SOURCE_DIR}/bin/runtest.sh ${CURTEST_DRIVER_ARGS})

      if (CURTEST_DEPENDS)
        set_tests_properties(${TestName} PROPERTIES DEPENDS "${CURTEST_DEPENDS}")
      endif()
      if (CURTEST_PROCESSORS)
        set_tests_properties(${TestName} PROPERTIES PROCESSORS "${CURTEST_PROCESSORS}")
      endif()
    endif()
  else ()
    if (NOT CURTEST_NO_COMPILE)
      add_test(${TestName} skip_test_dummy)
    endif()
  endif()
endmacro()
