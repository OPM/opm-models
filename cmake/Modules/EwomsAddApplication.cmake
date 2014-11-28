# add an application. the main difference to EwomsAddTest() is that
# applications are always compiled if CONDITION evaluates to true...
#
# EwomsAddApplication(AppName
#                     [EXE_NAME AppExecutableName]
#                     [CONDITION ConditionalExpression]
#                     [DRIVER_ARGS TestDriverScriptArguments]
#                     [SOURCES SourceFile1 SourceFile2 ...]
#                     [PROCESSORS NumberOfRequiredCoresForTest]
#                     [DEPENDS AppName1 AppName2 ...])
include(EwomsAddTest)

macro(EwomsAddApplication AppName)
  EwomsAddTest(${AppName} ${ARGN} ALWAYS_ENABLE)
endmacro()
