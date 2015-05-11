# add an application. the main difference to EwomsAddTest() is that
# applications are always compiled if CONDITION evaluates to true...
#
# EwomsAddApplication(AppName
#                     [EXE_NAME AppExecutableName]
#                     [CONDITION ConditionalExpression]
#                     [SOURCES SourceFile1 SourceFile2 ...])
include(EwomsAddTest)

macro(EwomsAddApplication AppName)
  EwomsAddTest(${AppName} ${ARGN} ALWAYS_ENABLE ONLY_COMPILE)
endmacro()
