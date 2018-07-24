#ifndef NONLINEAR_SYSTEM_H
#define NONLINEAR_SYSTEM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dlfcn.h>

#include "linear-system.h"
#include "problem.h"

// Flag to allow a function export ...
#define EXPORT_FN

// This works like a function pointer macro ...
#define SET_NONLINEAR_SYSTEM_SOLVER(name) EXPORT_FN double* name(double *J, double *f, const int n, set_linear_system_fn *solver_linear_system, set_problem_fn **functions)
typedef SET_NONLINEAR_SYSTEM_SOLVER(set_nonlinear_system_fn);

struct nonlinear_system_data
{
    void *handle;                       // Handle to the library that solves a linear system
    char *method_name;                  // Name of the solver method
    double *x;                          // Solution of the problem
    set_nonlinear_system_fn *solver;    // Pointer to the solver function
};

// Constructor and destructor
struct nonlinear_system_data* new_nonlinear_system (const int nonlinear_system_id);
void free_nonlinear_system (struct nonlinear_system_data *nls);

// Auxiliary functions 
set_nonlinear_system_fn* getNonlinearMethodFunction (void *handle, const int nonlinear_system_id);
char *getNonlinearMethodName (const int nonlinear_system_id);
//void printLinearSystem (struct linear_system_data *ls);
//void printMatrix (const char *name, const double *A, const int n);
//void printVector (const char *name, const double *v, const int n);

#endif