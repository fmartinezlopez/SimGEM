cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(TripleGemStack CXX)

# Default to an optimized build (-O3 -DNDEBUG). Without this, CMake uses no
# optimization flags and the avalanche loop runs several times slower.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Requires a Garfield++ install on CMAKE_PREFIX_PATH,
# e.g.  source <garfield_install>/share/Garfield/setupGarfield.sh
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

add_executable(generate_gas src/generate_gas.cpp)
target_link_libraries(generate_gas Garfield::Garfield)

add_executable(gain_stack src/gain_stack.cpp)
target_link_libraries(gain_stack Garfield::Garfield)

add_executable(vis_stack src/vis_stack.cpp)
target_link_libraries(vis_stack Garfield::Garfield)
