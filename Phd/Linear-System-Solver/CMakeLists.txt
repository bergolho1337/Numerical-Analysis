cmake_minimum_required(VERSION 2.8.9)
project(linear_system_solver)

# Bring the headers into the project
include_directories(include)

# The variable CMAKE_BINARY_DIR will be set to ${CMAKE_SOURCE_DIR}/bin
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR})

# EXECUTABLE_OUTPUT_PATH is the place where CMake will store the executables
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# This will include all the files on the /src folder that has the .cxx extension
file(GLOB SOURCES "src/*.cxx")

add_executable(linear_system_solver ${SOURCES}) 
