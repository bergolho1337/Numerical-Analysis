#include <iostream>
#include "matrix.h"
#include "conjGradient.h"

using namespace std;

int main ()
{
	// ENTRADA
	double *A;
	double *b;
	double *x;

	A = readMatrix_A();
	b = readVector_b();
	printMatrix(A);
	printVector(b);
	x = solveConjGradient(A,b);
	
}
