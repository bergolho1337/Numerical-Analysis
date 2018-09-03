#ifndef NONLINEAR_SYSTEM_SOLVER_H
#define NONLINEAR_SYSTEM_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Flag to allow a function export ...
#define EXPORT_FN

// This works like a function pointer macro ...
#define SET_LINEAR_SYSTEM_SOLVER(name) EXPORT_FN void name(double *A, double *b, const int n, double *x)
typedef SET_LINEAR_SYSTEM_SOLVER(set_linear_system_fn);

#define SET_PROBLEM_TYPE(name) EXPORT_FN double name(const double x[])
typedef SET_PROBLEM_TYPE(set_problem_fn);

const int MAX_ITER = 100;           // Maximum number of iterations
const double TOLERANCE = 1.0e-08;     // Tolerance of iterative method
const int REBUILD_JACOBIAN = 1;     // Rate which the Jacobian matrix will be rebuilted
const int GUESS_VECTOR = 1;         // Initial guess vector 
const double H = 0.0001;             // Size of the finite difference aproximation
const int FINITE_DIFFERENCE = 1;    // Type of finite difference

// Solver of Nonlinear Systems
double* NewtonFiniteDifference (double *J, double *f, const int n,\
                             set_linear_system_fn *solver_linear_system,\
                             set_problem_fn **functions);

// Auxiliary functions
void buildInitialGuess (double *x, const int n);
int hasConverged (const double *x, set_problem_fn **functions, const int n, const int iter);
int hasConverged2 (const double *s, const int n, const int iter);
double calcResidue (const double *x, const int n, set_problem_fn **functions);
double calcNorm (const double *v, const int n);
void buildJacobian_FiniteDifferences (double *J, double *x, const int n, set_problem_fn **functions);
void buildRHS (double *f, const double *x, set_problem_fn **functions, const int n);
void updateSolution (double *x, const double *s, const int n);
void printTypeFiniteDifference ();
void printCurrentSolution (const double *x, const int n, set_problem_fn **functions, const int iter);

#endif