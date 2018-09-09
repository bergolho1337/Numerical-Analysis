#ifndef LINEAR_SYSTEM_H
#define LINEAR_SYSTEM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dlfcn.h>

// Flag to allow a function export ...
#define EXPORT_FN

// This works like a function pointer macro ...
#define SET_LINEAR_SYSTEM_SOLVER(name) EXPORT_FN void name(double *A, double *b, const int n, double *x)
typedef SET_LINEAR_SYSTEM_SOLVER(set_linear_system_fn);

struct linear_system_data
{
    void *handle;                   // Handle to the library that solves a linear system
    char *method_name;              // Name of the solver method
    int n;                          // Size of the linear system
    double *A;                      // Coefficient matrix
    double *b;                      // Right-Hand-Side array
    double *x;                      // Solution array
    set_linear_system_fn *solver;   // Pointer to the solver function
};

// Constructor and destructor
struct linear_system_data* new_linear_system (const int linear_system_id);
void free_linear_system (struct linear_system_data *ls);

// Auxiliary functions 
set_linear_system_fn* getMethodFunction (void *handle, const int linear_system_id);
char *getMethodName (const int linear_system_id);
void printLinearSystem (struct linear_system_data *ls);
void printMatrix (const char *name, const double *A, const int n);
void printVector (const char *name, const double *v, const int n);

#endif