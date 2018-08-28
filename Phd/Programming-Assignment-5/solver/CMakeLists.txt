set(SOURCE_FILES_SOLVER solver.h solver.c
                        interpolation.h interpolation.c
                        newton_cotes.h newton_cotes.c)

ADD_LIBRARY(solver STATIC ${SOURCE_FILES_SOLVER})

set_property(TARGET solver PROPERTY POSITION_INDEPENDENT_CODE ON)