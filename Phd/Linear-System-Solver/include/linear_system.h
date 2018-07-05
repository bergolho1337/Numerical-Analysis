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

const int MAX_ITER = 500;           // Maximum number of iterations
const double EPSILON = 1.0e-05;     // Tolerance of the residue

class LinearSystem
{
private:
    unsigned int method;        // Method identifier
    unsigned int N;             // Size of the linear system
    double **A;                 // Pointer to the matrix
    double **LU;                // Pointer to the LU matrix
    double *b;                  // Pointer to the RHS
    double *x;                  // Pointer to the solution vector
    // Private functions
    void choosePivot (int &pivot_line, double &Amax, const int i);
    void reduceToRowEchelonForm ();
    void switchLines (int pivot[], const int pivot_line, const int i);
    void copyMatrixLU ();
public:
    LinearSystem (int argc, char *argv[]);
    ~LinearSystem ();
    void readMatrix (const char fname[]);
    double** allocMatrix (const int n, const int m);
    void printMatrix ();
    void readRHS (const char fname[]);
    void printRHS ();
    void allocRHS (const int n);
    double* allocVector (const int n);

    void solve ();

    void ForwardSubstitution ();
    void BackwardSubstitution ();
    void GaussianElimination ();
    void LUDecomposition ();

    void Jacobi ();
    void Gauss_Seidel ();
    void CG ();
    void BiCG ();
    void GMRES ();

    
    
    void checkSolution ();
    bool hasConverged (const int iter);
    double calcResidue ();

};

void Usage (const char pname[]);
void print (const char name[], double *v, const int N);

#endif