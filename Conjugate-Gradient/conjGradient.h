#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

double* solveConjGradient (double *A, double *b);
double calcLambda (double *A, double *r_old, double *aux);
double calcAlpha (double *r_old, double *r);
void calcP (double *r, double alpha, double *p_old, double *aux);
double calcBeta (double *r, double *A, double *p, double *aux);
void calcV (double *v, double *v_old, double beta, double *p, double *aux);
void calcR (double *r, double *r_old, double beta, double *A, double *p, double *aux);
bool checkConvergence (double *r, double *r_old, double *aux, int N);
