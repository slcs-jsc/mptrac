# Minimum version of CMake
cmake_minimum_required(VERSION 3.10)

#Project name and version
project(MPTRAC_DD VERSION 0.0)

# Collect all .c files from the src directory
file(GLOB SOURCES "src/*.c")

# Add the library
add_library(MPTRAC_DD STATIC ${SOURCES})

# Specifcy include directories
target_include_directories(MPTRAC_DD PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>$<INSTALL_INTERFACE:src>)
