#ifndef LINEAR_SYSTEM_SOLVER_H
#define LINEAR_SYSTEM_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int MAX_ITER = 500;           // Maximum number of iterations
const double EPSILON = 1.0e-08;     // Tolerance of iterative method
const int PRECONDITIONER = 0;       // Flag to set the Jacobi preconditioner

// Solver of Linear Systems
void Jacobi (double *A, double *b, const int n, double *x);
void Gauss_Seidel (double *A, double *b, const int n, double *x);
void CG (double *A, double *b, const int n, double *x);
void BiCG (double *A, double *b, const int n, double *x);
void LUDecomposition (double *A, double *b, const int n, double *x);

// Auxiliary functions
void swap (double **a, double **b);
int hasConverged (const double *A, const double *b, const double *x, const int n, const int iter);
double calcResidue (const double *A, const double *b, const double *x, const int n);
void checkSolution (const double *A, const double *b, const double *x, const int n, const int iter);
void choosePivot (double *LU, int *pivot_line, double *Amax, const int i, const int N);
void switchLines (double *LU, int pivot[], const int pivot_line, const int i, const int N);

#endif