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
    void *handle;
    char *method_name;
    int n;
    double *A;
    double *b;
    double *x;
    set_linear_system_fn *solver;
};

struct linear_system_data* new_linear_system (const int linear_system_id);
void free_linear_system (struct linear_system_data *ls);

set_linear_system_fn* getMethodFunction (void *handle, const int linear_system_id);
char *getMethodName (const int linear_system_id);
void printLinearSystem (struct linear_system_data *ls);
void printMatrix (const char *name, const double *A, const int n);
void printVector (const char *name, const double *v, const int n);

#endif