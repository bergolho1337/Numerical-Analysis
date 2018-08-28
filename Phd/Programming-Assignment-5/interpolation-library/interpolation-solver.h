#ifndef INTERPOLATION_SOLVER_H
#define INTERPOLATION_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int NEVAL = 100;              // Number of evaluations

// Interpolation polynomials
double Lagrange (double *x, double *y, const int a, const int b, const double z);
double Newton (double *x, double *y, const int a, const int b, const double z);

// Auxiliary functions


#endif