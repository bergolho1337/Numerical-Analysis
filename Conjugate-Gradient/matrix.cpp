#include "matrix.h"

// Size of the system
int N;			

// Sum two vectors 'a' and 'b' and stores in 'c'
void sumVectorVector (double *a, double *b, double *c)
{
	int i;
	for (i = 0; i < N; i++)
		c[i] = a[i]+b[i];
}

// Multiply the vector 'a' by a scalar 's'
void multiplyScalarVector (double *a, double s)
{
	int i;
	for (i = 0; i < N; i++)
		a[i] *= s;

}

double multiplyVectorVector (double *a, double *b)
{
	int i;
	double result = 0.0;
	for (i = 0; i < N; i++)
		result += a[i]*b[i];
	return (result);
}

// Multiply the matrix 'A' and the vector 'b' and stores in the vector 'c' 
void multiplyMatrixVector (double *A, double *b, double *c)
{
	int i, j;
	double value;
	for (i = 0; i < N; i++)
	{
		value = 0.0;
		for (j = 0; j < N; j++)
			value += A[i*N+j]*b[j];
		c[i] = value;
	}
}

double* readMatrix_A ()
{
	double *A;
	int i, j;
	
	cin >> N;
	// Save the order of the matrix
	A = new double[N*N];
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			cin >> A[i*N+j];
	return (A);
	
}

double* readVector_b ()
{
	double *b;
	int i;

	b = new double[N];
	for (i = 0; i < N; i++)
		cin >> b[i];
	return (b);
}

void printMatrix (double *A)
{
	int i, j;
	cout << "[!] Printing matrix A" << endl;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
			cout << A[i*N+j] << " ";
		cout << endl;
	}
	cout << endl;
}

void printVector (double *b)
{
	int i;
	cout << "[!] Printing vector b" << endl;
	for (i = 0; i < N; i++)
		cout << b[i] << endl;
	cout << endl;
}

void copyVector (double *a, double *b)
{
	int i;
	for (i = 0; i < N; i++)
		b[i] = a[i];
}

int getSize ()
{
	return (N);
}

