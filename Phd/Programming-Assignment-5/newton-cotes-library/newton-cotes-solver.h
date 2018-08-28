#ifndef NEWTON_COTES_SOLVER_H
#define NEWTON_COTES_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../solver/newton_cotes.h"

// Number of subintervals (Needs to be multiple of 6 to work with all the 3 methods ...)
const int NSUBINTERVAL = 3600;              

// Newton-Cotes rules
double Trapezium (const int minid, const int maxid, struct interpolation_data *idata, const int degree);
double Simpson13 (const int minid, const int maxid, struct interpolation_data *idata, const int degree);
double Simpson38 (const int minid, const int maxid, struct interpolation_data *idata, const int degree);
double Mixed (const int minid, const int maxid, struct interpolation_data *idata, const int degree);

// Auxiliary functions


#endif