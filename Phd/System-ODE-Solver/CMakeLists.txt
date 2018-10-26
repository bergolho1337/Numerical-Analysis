cmake_minimum_required(VERSION 3.0.1)
project(solver)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

ADD_EXECUTABLE(solver main.cpp)

TARGET_LINK_LIBRARIES(solver m)
