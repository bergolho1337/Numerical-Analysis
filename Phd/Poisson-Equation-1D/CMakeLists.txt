cmake_minimum_required(VERSION 3.0.1)
project(PoissonEquation1d)

IF (CMAKE_VERSION VERSION_LESS "3.1")
    IF (CMAKE_C_COMPILER_ID STREQUAL "GNU")
        SET(CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
    ENDIF(CMAKE_C_COMPILER_ID)
ELSE()
    SET(CMAKE_C_STANDARD 99)
ENDIF()

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-strict-aliasing -Wall -Wno-unused-function")

#ADD_SUBDIRECTORY(src/utils)
#ADD_SUBDIRECTORY(src/sparse_matrix)
#ADD_SUBDIRECTORY(src/config)

ADD_EXECUTABLE(PoissonEquation1d src/main.cpp)

#ADD_DEPENDENCIES(PoissonEquation1d utils sparse_matrix config )

#TARGET_LINK_LIBRARIES(PoissonEquation1d utils sparse_matrix config dl m)

###