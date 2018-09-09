SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT LINEAR SYSTEM SOLVER LIBRARY ##############################
ADD_LIBRARY(default-least-squares-solver ${LIB_TYPE} least-squares-solver.c)
TARGET_LINK_LIBRARIES(default-least-squares-solver)
##################################################################################
