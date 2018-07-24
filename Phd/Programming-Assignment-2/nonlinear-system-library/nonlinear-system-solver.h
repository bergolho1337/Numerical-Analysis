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

const int MAX_ITER = 3;           // Maximum number of iterations
const double EPSILON = 1.0e-05;     // Tolerance of iterative method
const int REBUILD_JACOBIAN = 1;     // Rate which the Jacobian matrix will be rebuilted
const int GUESS_VECTOR = 2;         // Initial guess vector 
const double H = 0.001;             // Size of the finite difference aproximation
const int FINITE_DIFFERENCE = 1;     // Type of finite difference

// Solver of Nonlinear Systems
void NewtonFiniteDifference (double *J, double *f, const int n, double *x,\
                             set_linear_system_fn *solver_linear_system,\
                             set_problem_fn **functions);
//void Gauss_Seidel (double *A, double *b, const int n, double *x);
//void CG (double *A, double *b, const int n, double *x);
//void BiCG (double *A, double *b, const int n, double *x);

// Auxiliary functions
void buildInitialGuess (double *x, const int n);
int hasConverged (const double *x, set_problem_fn **functions, const int n, const int iter);
double calcResidue (const double *x, const int n, set_problem_fn **functions);
void buildJacobian_FiniteDifferences (double *J, double *x, const int n, set_problem_fn **functions);
void printTypeFiniteDifference ();
//void swap (double **a, double **b);
//int hasConverged (const double *A, const double *b, const double *x, const int n, const int iter);
//double calcResidue (const double *A, const double *b, const double *x, const int n);
//void checkSolution (const double *A, const double *b, const double *x, const int n, const int iter);

#endif