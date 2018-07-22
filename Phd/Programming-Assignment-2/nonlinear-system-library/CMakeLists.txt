SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT LINEAR SYSTEM SOLVER LIBRARY ##############################
SET(LINEAR_SYSTEM_SOLVER_SRC ../linear-system-library/linear-system-solver.c ../linear-system-library/linear-system-solver.h)
ADD_LIBRARY(default-nonlinear-system-solver ${LIB_TYPE} ${LINEAR_SYSTEM_SOLVER_SRC} nonlinear-system-solver.c)
TARGET_LINK_LIBRARIES(default-nonlinear-system-solver)
##################################################################################
