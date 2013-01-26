# set the CMake module include path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# eWoms unit-testing infrastructure
include(EwomsAddTest)

EwomsAddTest(test_adaptive2p2c CONDITION ALUGRID_FOUND)
EwomsAddTest(test_dec2p2c)
EwomsAddTest(test_multiphysics2p2c 
  DRIVER_ARGS --simulation test_multiphysics2p2c --end-time=2000)
