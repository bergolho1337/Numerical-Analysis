#ifndef LINEAR_SYSTEM_SOLVER_H
#define LINEAR_SYSTEM_SOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int MAX_ITER = 500;
const double EPSILON = 1.0e-05;
const int PRECONDITIONER = 0; 

void Jacobi (double *A, double *b, const int n, double *x);
void Gauss_Seidel (double *A, double *b, const int n, double *x);
void CG (double *A, double *b, const int n, double *x);
void BiCG (double *A, double *b, const int n, double *x);

void swap (double **a, double **b);
int hasConverged (const double *A, const double *b, const double *x, const int n, const int iter);
double calcResidue (const double *A, const double *b, const double *x, const int n);
void checkSolution (const double *A, const double *b, const double *x, const int n, const int iter);

#endif