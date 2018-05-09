#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* setMatrix ()
{
	double *A = (double*)malloc(sizeof(double)*3*3);
	A[0] = 4;
	A[1] = -6;
	A[2] = 2;
	A[3] = -6;
	A[4] = 10;
	A[5] = -5;
	A[6] = 2;
	A[7] = -5;
	A[8] = 30;
	return A;
}

double* setVector ()
{
	double *b = (double*)calloc(3,sizeof(double));
	return b;
}

void printMatrix (double *A)
{
	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
			printf("%e ",A[i*3+j]);
		printf("\n");
	}
	printf("\n");
}

void Cholesky (double *A)
{
	int i, j, k;
	int condErro = 0;
	double soma, t, r;
	for (j = 0; j < 3; j++)
	{
		soma = 0;
		for (k = 0; k <= j-1; k++)
		{
			soma += pow(A[j*3+k],2);
		}
		t = A[j*3+j] - soma;
		if (t > 0)
		{
			A[j*3+j] = sqrt(t);
			r = 1 / A[j*3+j];
		}
		else
		{
			condErro = 1;
			printf("[-] A matriz nao eh positiva definida.\n");
			exit(1);
		}
		for (i = j+1; i < 3; i++)
		{
			soma = 0;
			for (k = 0; k <= j-1; k++)
			{
				soma += A[i*3+k]*A[j*3+k];
			}
			A[i*3+j] = (A[i*3+j] - soma)*r;
			A[j*3+i] = A[i*3+j];
		}
	}
}

void solveLinearSystem (double *A, double *b, double *x)
{
	int i, j, k;
	double *y = (double*)malloc(sizeof(double)*3);
	double soma;
	y[0] = b[0] / A[0];
	for (i = 1; i < 3; i++)
	{
		soma = 0;
		for (j = 0; j < i; j++)
		{
			soma += A[i*3+j]*y[j];
		}
		y[i] = (b[i] - soma) / A[i*3+i];
	}
	printf("Vetor y:\n");
	for (i = 0; i < 3; i++)
		printf("%lf\n",y[i]);
	x[2] = y[2] / A[2*3+2];
	for (i = 1; i >= 0; i--)
	{
		soma = 0;
		for (j = 2; j > i; j--)
		{
			soma += A[i*3+j]*x[j];
		}
		x[i] = (y[i] - soma) / A[i*3+i];
	}
	printf("\nSolucao:\n");
	for (i = 0; i < 3; i++)
		printf("%lf\n",x[i]);
}

int main ()
{
	double *A, *b, *x;
	// Alocar estruturas
	A = setMatrix();
	b = setVector();
	x = setVector();

	b[0] = 1;
	Cholesky(A);
	printMatrix(A);
	solveLinearSystem(A,b,x);
	free(A);
	return 0;
}
