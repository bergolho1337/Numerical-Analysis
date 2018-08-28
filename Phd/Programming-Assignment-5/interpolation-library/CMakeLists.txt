SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT INTERPOLATION LIBRARY #####################################
ADD_LIBRARY(default-interpolation-solver ${LIB_TYPE} interpolation-solver.c)
TARGET_LINK_LIBRARIES(default-interpolation-solver)
##################################################################################