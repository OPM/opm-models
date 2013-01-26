# set the CMake module include path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# eWoms unit-testing infrastructure
include(EwomsAddTest)

EwomsAddTest(test_propertysystem DRIVER_ARGS --plain test_propertysystem)
