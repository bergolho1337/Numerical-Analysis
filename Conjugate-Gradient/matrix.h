#include <iostream>
#include <cstdlib>

using namespace std;

double* readMatrix_A ();
double* readVector_b ();
void printMatrix (double *A);
void printVector (double *b);
void sumVectorVector (double *a, double *b, double *c);
void multiplyMatrixVector (double *A, double *b, double *c);
double multiplyVectorVector (double *a, double *b);
void multiplyScalarVector (double *a, double s);
void copyVector (double *a, double *b);
int getSize ();
