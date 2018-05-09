#include "matriz.h"

void imprimeMatriz (double **m, int n)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		cout << endl;
		for (j = 0; j < n; j++)
			printf("%e\t", m[i][j]);
	}
	cout << endl;
}

void imprimeVetor (double *v, int n)
{
	int i;
	for (i = 0; i < n; i++)
		cout << v[i] << endl;
}

double** multiplicaMatrizMatriz (double **a, double **b, int n)
{
	int i, j, k;
	double **c = new double*[n];
	double soma;
	for (i = 0; i < n; i++)
		c[i] = new double[n];
	for (i = 0; i < n; i++)
	{	
		for (j = 0; j < n; j++)
		{
			soma = 0;			
			for (k = 0; k < 6; k++)			
				soma = soma + (a[i][k]*b[k][j]);
			c[i][j] = soma;
		}
	}
	return (c);		
}

double* multiplicaMatrizVetor (double **a, double *b, int n)
{
	int i, j;
	double *c = new double[n];
	double soma;
	for (i = 0; i < n; i++)
	{	
		soma = 0;		
		for (j = 0; j < n; j++)
		{						
			soma = soma + (a[i][j]*b[j]);
		}
		c[i] = soma;
	}
	return (c);		
}

double* somaVetorVetor (double *a, double *b, int n)
{
	int i;
	double *c = new double[n];
	for (i = 0; i < n; i++)
	{	
		c[i] = a[i]+b[i];
	}
	return (c);
}

double multiplicaVetorVetor (double *a, double *b, int n)
{
	int i;
	double soma = 0;
	for (i = 0; i < n; i++)
		soma += (a[i]*b[i]);
	return (soma);
}

void multiplicaEscalarVetor (double *a, double alpha, int n)
{
	int i;
	for (i = 0; i < n; i++)
		a[i] = a[i]*alpha;
}

// Vetor 'a' = fonte e vetor 'b' = destino
void copiaVetor (double *a, double *b, int n)
{
	int i;
	for (i = 0; i < n; i++)
		b[i] = a[i];
}
