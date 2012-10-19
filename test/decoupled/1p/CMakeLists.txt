# add build targets
add_executable("test_diffusion" test_diffusion.cc)
target_link_libraries("test_diffusion" ${EwomsLinkLibraries})

add_executable("test_dec1p" test_1p.cc)
target_link_libraries("test_dec1p" ${EwomsLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${EwomsLinkDirectories})
include_directories(${EwomsIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_diffusion"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

