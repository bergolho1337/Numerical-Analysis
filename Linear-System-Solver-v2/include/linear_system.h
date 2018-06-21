#ifndef __LINEAR_SYSTEM_H__
#define __LINEAR_SYSTEM_H__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>

#define PRINT_LINE "=================================================================================================================="

const int MAX_ITER = 500;           // Maximum number of iterations
const double EPSILON = 1.0e-05;     // Tolerance of the residue

class LinearSystem
{
private:
    unsigned int method;        // Method identifier
    unsigned int N;             // Size of the linear system
    double **A;                 // Pointer to the matrix
    double *b;                  // Pointer to the RHS
    double *x;                  // Pointer to the solution vector
public:
    LinearSystem (int argc, char *argv[]);
    ~LinearSystem ();
    void readMatrix (const char fname[]);
    void allocMatrix (const int n, const int m);
    void printMatrix ();
    void readRHS (const char fname[]);
    void printRHS ();
    void allocRHS (const int n);

    void solve ();
    void Jacobi ();
    void Gauss_Seidel ();
    void CG ();
    void BCG ();

    void checkSolution ();
    bool hasConverged (const int iter);
    double calcResidue ();

};

void Usage (const char pname[]);
void print (const char name[], double *v, const int N);

#endif