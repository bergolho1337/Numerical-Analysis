#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <iostream>

#include "lagrange.h"
#include "newton.h"
#include "csplines.h"

using namespace std;

// Number of evaluations of the polynomial
const int NEVAL = 100;

// Flag to allow a function export ...
#define EXPORT_FN

// This works like a function pointer macro ...
#define SET_INTERPOLATION(name) EXPORT_FN double name(double *x, double *y,\
                                                    const int a, const int b, const int degree,\
                                                    const double z)
typedef SET_INTERPOLATION(set_interpolation_fn);

set_interpolation_fn* get_interpolation_function (const int id);

#endif



