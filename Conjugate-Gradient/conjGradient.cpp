#include "conjGradient.h"
#include "matrix.h"

double Toler = 1.0e-02;

double* solveConjGradient (double *A, double *b)
{
	int i, N, Iter;

	N = getSize();
	Iter = 2;
	double *v_old = new double[N];
	double *v = new double[N];
	double *p_old = new double[N];
	double *p = new double[N];
	double *r_old = new double[N];
	double *r = new double[N];
	double *aux = new double[N];
	double lambda, alpha, beta;

	// Initialize the solutionwith zeros
	for (i = 0; i < N; i++)
		v_old[i] = 0.0;

	multiplyMatrixVector(A,v_old,aux);
	sumVectorVector(aux,b,r_old);
	copyVector(r_old,p_old);
	multiplyScalarVector(r_old,-1);
	lambda = calcLambda(A,r_old,aux);

	copyVector(p,aux);
	multiplyScalarVector(aux,lambda);
	sumVectorVector(v_old,aux,v);
	multiplyMatrixVector(A,p_old,aux);
	multiplyScalarVector(aux,lambda);
	sumVectorVector(r_old,aux,r);

	do
	{
		alpha = calcAlpha(r_old,r);
		calcP(r,alpha,p_old,aux);
		beta = calcBeta(r,A,p_old,aux);
		calcV(v,v_old,beta,p_old,aux);
		calcR(r,r_old,beta,A,p_old,aux);
		copyVector(v,v_old);
		copyVector(r,r_old);
		Iter++;	
	} while (Iter < 4);

	delete [] v_old;
	delete [] p;
	delete [] p_old;
	delete [] r_old;
	delete [] r;

	cout << "Solution" << endl;
	for (i = 0; i < N; i++)
		cout << v[i] << endl;
	cout << "Number of iterations = " << Iter << endl;
	return (v);
}

bool checkConvergence (double *r, double *r_old, double *aux, int N)
{
	int i;
	double max_value_Num = 0.0;
	double max_value_Den = 0.0;
	for (i = 0; i < N; i++)
	{
		aux[i] = r[i]-r_old[i];
		if (aux[i] > max_value_Num)
			max_value_Num = aux[i];
		if (r[i] > max_value_Den)
			max_value_Den = r[i];
	}
	if (max_value_Num/max_value_Den < Toler)
		return (true);
	else
		return (false);
	
}

double calcLambda (double *A, double *r_old, double *aux)
{
	double Num, Den;
	Num = multiplyVectorVector(r_old,r_old);
	multiplyMatrixVector(A,r_old,aux);
	Den = multiplyVectorVector(aux,r_old);
	return (Num/Den);
}

double calcAlpha (double *r_old, double *r)
{
	double Num, Den;
	Num = multiplyVectorVector(r,r);
	Den = multiplyVectorVector(r_old,r_old);
	return (Num/Den);
}

void calcP (double *r, double alpha, double *p_old, double *aux)
{
	copyVector(r,aux);
	multiplyScalarVector(aux,-1);
	multiplyScalarVector(p_old,alpha);
	sumVectorVector(aux,p_old,p_old);
}

double calcBeta (double *r, double *A, double *p, double *aux)
{
	double Num, Den;
	Num = multiplyVectorVector(r,r);
	multiplyMatrixVector(A,p,aux);
	Den = multiplyVectorVector(aux,p);
	return (Num/Den);
}

void calcV (double *v, double *v_old, double beta, double *p, double *aux)
{
	copyVector(p,aux);
	multiplyScalarVector(aux,beta);
	sumVectorVector(v_old,aux,v);
}

void calcR (double *r, double *r_old, double beta, double *A, double *p, double *aux)
{
	multiplyMatrixVector(A,p,aux);
	multiplyScalarVector(aux,beta);
	sumVectorVector(r_old,aux,r);
}
