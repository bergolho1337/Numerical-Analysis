SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT NEWTON-COTES LIBRARY ######################################
ADD_LIBRARY(default-newton-cotes-solver ${LIB_TYPE} newton-cotes-solver.c)
TARGET_LINK_LIBRARIES(default-newton-cotes-solver)
##################################################################################