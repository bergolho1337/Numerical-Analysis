SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT NONLINEAR SYSTEM SOLVER LIBRARY ##############################
SET(NONLINEAR_SYSTEM_SOLVER_SRC ../nonlinear-system-library/nonlinear-system-solver.c ../nonlinear-system-library/nonlinear-system-solver.h)
ADD_LIBRARY(default-nonlinear-system-solver ${LIB_TYPE} ${NONLINEAR_SYSTEM_SOLVER_SRC} nonlinear-system-solver.c)
TARGET_LINK_LIBRARIES(default-nonlinear-system-solver)
##################################################################################
