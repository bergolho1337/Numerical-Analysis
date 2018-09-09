SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT LINEAR SYSTEM SOLVER LIBRARY ##############################
ADD_LIBRARY(default-linear-system-solver ${LIB_TYPE} linear-system-solver.c)
TARGET_LINK_LIBRARIES(default-linear-system-solver)
##################################################################################
