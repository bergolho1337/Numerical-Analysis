set(SOURCE_FILES_SOLVER solver.h solver.c
                        problem.h problem.c)

ADD_LIBRARY(solver STATIC ${SOURCE_FILES_SOLVER})

set_property(TARGET solver PROPERTY POSITION_INDEPENDENT_CODE ON)