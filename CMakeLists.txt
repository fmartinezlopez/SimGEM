cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(SimGEM CXX)

# Default to an optimized build (-O3 -DNDEBUG). Without this, CMake uses no
# optimization flags and the avalanche loop runs several times slower.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Requires a Garfield++ install on CMAKE_PREFIX_PATH,
# e.g.  source <garfield_install>/share/Garfield/setupGarfield.sh
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

foreach(prog gain vis transport gastable penning counter)
  add_executable(simgem_${prog} src/${prog}.cpp)
  target_include_directories(simgem_${prog} PRIVATE src)
  target_link_libraries(simgem_${prog} Garfield::Garfield)
endforeach()
