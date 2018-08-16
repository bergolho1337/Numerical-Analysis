#ifndef __LINEAR_SYSTEM_H__
#define __LINEAR_SYSTEM_H__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstring>

#define PRINT_LINE "=================================================================================================================="
#define PRINT_DANGER "****************************************************************************************************"

const int MAX_ITER = 500;           // Maximum number of iterations
const double EPSILON = 1.0e-08;     // Tolerance of the residue

class LinearSystem
{
private:
    unsigned int method;        // Method identifier
    unsigned int N;             // Size of the linear system
    double **A;                 // Pointer to the matrix
    double **LU;                // Pointer to the LU matrix
    double *b;                  // Pointer to the RHS
    double *x;                  // Pointer to the solution array
    double *r;                  // Pointer to the residue array
    double *d;                  // Pointer to the correction array
public:
    LinearSystem (int argc, char *argv[]);
    ~LinearSystem ();
    void readMatrix (const char fname[]);
    double* allocVector (const int n);
    void printMatrix ();
    void readRHS (const char fname[]);
    void printRHS ();

    void solve ();              // Solve the linear system using the specified method
    void cond ();               // Calculate the condition number of the system
    void write ();              // Write the solution of the system to a file

};

void ForwardSubstitution (double **A, double *b, double *x, const int N);
void BackwardSubstitution (double **A, double *b, double *x, const int N);
void GaussianElimination (double **A, double *b, double *x, const int N);
void LUDecomposition (double **A, double **LU, double *b, double *x, const int N);

void Jacobi (double **A, double *b, double *x, const int N);
void Gauss_Seidel (double **A, double *b, double *x, const int N);
void CG (double **A, double *b, double *x, const int N);
void BiCG (double **A, double *b, double *x, const int N);
void GMRES (double **A, double *b, double *x, const int N);

bool checkSolution (double **A, const double *x, const double *b, const int N);
bool refineSolution (double **A, double *x, double *b, const int N);
bool hasConverged (double **A, const double *x, const double *b, const int iter, const int N);
double calcResidue (double **A, const double *x, const double *b, const int N);
double calcResidueInf (double **A, const double *x, const double *b, const int N);
void compResidue (double **A, const double *x, const double *b, double *r, const int N);

void choosePivot (double **LU, int &pivot_line, double &Amax, const int i, const int N);
void reduceToRowEchelonForm (double **A, double *b, const int N);
void switchLines (double **LU, int pivot[], const int pivot_line, const int i, const int N);
void copyMatrixLU (double **A, double **LU, const int N);
void buildElementaryVector (double *e, const int k, const int N);
void insertColumnMatrix (double **A, double *v, const int k, const int N);

void Usage (const char pname[]);
void print (const char name[], double *v, const int N);

double** allocMatrix (const int n, const int m);

void MatrixMatrixMultiplication (double **A, double **B, const int N);
double calcMatrixNormInf (double **A, const int N);

#endif