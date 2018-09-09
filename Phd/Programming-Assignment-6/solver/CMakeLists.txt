set(SOURCE_FILES_SOLVER solver.h solver.c
                        least-squares.h least-squares.c
                        linear-system.h linear-system.c)

ADD_LIBRARY(solver STATIC ${SOURCE_FILES_SOLVER})

set_property(TARGET solver PROPERTY POSITION_INDEPENDENT_CODE ON)