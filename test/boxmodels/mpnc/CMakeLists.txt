# build target for the twophase-twocomponent test problem
add_executable("test_mpnc" test_mpnc.cc)

add_executable("obstacle_2p2c" test_mpnc.cc)
set_target_properties("obstacle_2p2c" PROPERTIES COMPILE_FLAGS "-DUSE_2P2C=1")

# add required libraries and includes to the build flags 
link_directories(${DumuxLinkDirectories})
include_directories(${DumuxIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_mpnc"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

