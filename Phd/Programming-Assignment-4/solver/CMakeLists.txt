set(SOURCE_FILES_SOLVER solver.h solver.cpp
                        lagrange.h lagrange.cpp
                        newton.h newton.cpp
                        csplines.h csplines.cpp)

ADD_LIBRARY(solver STATIC ${SOURCE_FILES_SOLVER})

set_property(TARGET solver PROPERTY POSITION_INDEPENDENT_CODE ON)