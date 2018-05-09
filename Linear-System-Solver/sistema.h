#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "matriz.h"

using namespace std;

double** leMatriz ();
double* leVetor ();
void imprimeSolucao (double *x);

// Resolve um sistema linear pela substituição LU.
double* LU (double **A, double *b);

// Resolve um sistema linear por Cholesky, a matriz [A] deve ser simétrica e positiva definida.
double* Cholesky (double **A, double *b);

// Checa se a matriz dos coeficientes é positiva definida
bool ehPositivaDefinida (double **A);

// Resolve um sistema pelo método iterativo de Gauss-Seidel.
double* Gauss_Seidel (double **A, double *b);
bool testeConvergencia (double **A);
void imprimeIteracao (int Iter, double *x, double NormaRel);

// Resolve um sistema pelo método dos Gradientes Conjugados
double* GradientesConjugados (double **A, double *b);
double calcLambda (double **A, double *r, double *d);
double* calcX (double *x_ant, double lambda, double *d_ant);
double* calcR (double **A, double *r_ant, double *d, double lambda);
double calcBeta (double *r, double *r_ant);
double* calcD (double *r, double beta, double *d_ant);


