# set the CMake module include path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# eWoms unit-testing infrastructure
include(EwomsAddTest)

EwomsAddTest(test_diffusion DRIVER_ARGS --simulation-diffusion test_diffusion 3)
EwomsAddTest(test_1p DRIVER_ARGS --plain test_1p)
