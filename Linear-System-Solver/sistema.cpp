#include "sistema.h"

int n;					// Ordem da matriz.
FILE *arq = fopen("entrada4.txt","r");

// Resolve um sistema linear pela substituição LU.
double* LU (double **A, double *b)
{
	int i, j, k, p;
	double *pivot = new double[n];
	double Amax, t, m, r, Mult;
	// 1 PASSO: Transformar a matriz A do problema em duas matrizes triangulares L e U. 	
	for (i = 0; i < n; i++)
		pivot[i] = i;
	for (j = 0; j < n-1; j++)
	{
		cout << "Iteracao " << j << endl;
		imprimeMatriz(A,n);		
		// Escolher pivot
		p = j;
		Amax = abs(A[j][j]);
		// Verifica na coluna a ser eliminada qual elemento possui o maior valor absoluto, este elemento será o pivô.		
		for (k = j+1; k < n; k++)
		{
			if (abs(A[k][j]) > Amax)
			{
				Amax = abs(A[k][j]);
				p = k;
			}
		}
		// Se (p != j) então deve-se trocar de linhas 
		if (p != j)
		{
			for (k = 0; k < n; k++)
			{
				t = A[j][k];
				A[j][k] = A[p][k];
				A[p][k] = t;	
			}
			m = pivot[j];
			pivot[j] = pivot[p];
			pivot[p] = m;
		}
		if (abs(A[j][j]) != 0)
		{
			// Eliminação de Gauss
			r = 1 / A[j][j];
			for (i = j+1; i < n; i++)
			{
				Mult = A[i][j]*r;
				A[i][j] = Mult;
				for (k = j+1; k < n; k++)
					A[i][k] = A[i][k] - Mult*A[j][k];
			}
		}
	}
	cout << "\nDecomposicao L/U" << endl;
	imprimeMatriz(A,n);
	// A matriz A agora é L\U. Elementos abaixo da diagonal principal são de L, os da diagonal principal para cima pertencem a U.
	
	// 2 PASSO: Realizar substituições sucessivas para resolver o sistema triangular inferior: Ly = b		
	double *y = new double[n];
	double soma;
	k = pivot[0];
	y[0] = b[k];
	// Percorre somente os elementos de L em A.
	for (i = 1; i < n; i++)
	{
		soma = 0;
		for (j = 0; j <= i-1; j++)
		{
			soma += A[i][j]*y[j]; 
		}
		k = pivot[i];
		y[i] = b[k] - soma;	
	}
	// 3 PASSO: Realizar substituições retroativas para resolver o sistema triangular superior: Ux = y 			
	double *x = new double[n];
	x[n-1] = y[n-1] / A[n-1][n-1];
	for (i = n-2; i >= 0; i--)
	{
		soma = 0;
		for (j = i+1; j < n; j++)
			soma += A[i][j]*x[j];
		x[i] = (y[i] - soma) / A[i][i];
	}	

	delete [] pivot;
	delete [] y;

	imprimeSolucao(x);
	return (x);
}

bool ehPositivaDefinida (double **A)
{
	int i, j;
	double soma, t;
	bool tick = true;	
	for (i = 0; i < n; i++)
	{
		soma = 0;
		for (j = 0; j < i; j++)
			soma += (A[i][j]*A[i][j]);
		t = A[i][i] - soma;
		if (t <= 0)			// Não é positiva definida
			tick = false;
	}
	return (tick);
}

double* Cholesky (double **A, double *b)
{
	int i, j, k;
	bool condErro = false;
	double Det = 1, soma, t, r;
	for (i = 0; i < n; i++)
	{
		soma = 0;
		for (j = 0; j < i; j++)
			soma = soma + (A[i][j]*A[i][j]);
		t = A[i][i] - soma;
		if (t > 0)			// É positiva definida
		{
			A[i][i] = sqrt(t);	// Elemento da diagonal principal
			r = 1/A[i][i];
			Det = Det*t;
		}	
		else
		{
			condErro = true;
			cout << "[-] ERRO! A matriz [A] nao eh definida positiva" << endl;
			exit(1);
		}
		for (k = i+1; k < n; k++)
		{
			soma = 0;
			for (j = 0; j < i; j++)
				soma = soma + (A[k][j]*A[i][j]);
			A[k][i] = (A[k][i] - soma)*r;
		}
	}
	for (i = 0; i < n; i++)
	{
		for (j = i+1; j < n; j++)
			A[i][j] = A[j][i];
	}
	// Matriz [A] agora é simétrica.	
	
	double *y = new double[n];
	double *x = new double[n];
	// Substituições sucessivas
	for (i = 0; i < n; i++)
	{
		soma = 0;		
		for (j = 0; j < i; j++)
		{
			soma += A[i][j]*y[j];
		}
		y[i] = (b[i] - soma)/A[i][i];
	}
	// Substituições retroativas
	for (i = n-1; i >= 0; i--)
	{
		soma = 0;		
		for (j = i; j < n; j++)
		{
			soma += A[j][i]*x[j];
		}
		x[i] = (y[i] - soma)/A[i][i];
	}
	imprimeSolucao(x);
	return (x);
}

void imprimeIteracao (int Iter, double *x, double NormaRel)
{
	int i;
	cout << "\n=== ITERACAO " << Iter << " ===" << endl;
	cout << "Norma relativa = " << NormaRel << endl;
	cout << "--> x" << endl;	
	for (i = 0; i < n; i++)
		cout << x[i] << endl;
}

// Testa se a matriz [A] é diagonal estritamente dominante
bool testeConvergencia (double **A)
{
	int i, j;
	double soma;
	for (i = 0; i < n; i++)
	{
		soma = 0;		
		for (j = 0; j < n; j++)		
		{
			if (i != j)
				soma += abs(A[i][j]);
		}
		if (abs(A[i][i]) < soma)
			return (false);
	}
	return (true);	
}

double* Gauss_Seidel (double **A, double *b)
{
	int i, j;
	double NormaNum, NormaDen;		// Numerador e denominador para o cálculo do erro relativo (< Tolerancia)
	double Toler, IterMax, r, Iter, NormaRel, soma, t;
	double *x = new double[n];
	double *v = new double[n]; 		// Vetor auxiliar para armazenar os passos anteriores
	Toler = 1.0e-05;
	IterMax = 50;
	// Checar critério de convergência
	if (!testeConvergencia(A))
	{
		cout << "[-]ERRO! A matriz [A] nao passou no teste de convergencia!" << endl;
		exit(1);
	}
	// Construção das matrizes de iteração
	for (i = 0; i < n; i++)
	{
		r = 1/A[i][i];
		for (j = 0; j < n; j++)
		{
			if (i != j)
				A[i][j] = A[i][j]*r;
		}
		b[i] = b[i]*r;
		x[i] = b[i];		
	}
	Iter = 0;
	// Iterações de Gauss-Seidel
	while ( ( NormaRel >= Toler ) || ( Iter <= IterMax ) )		// Teste de convergência
	{
		Iter++;
		for (i = 0; i < n; i++)
		{
			soma = 0;
			for (j = 0; j < n; j++)
			{
				if (i != j)
					soma = soma + (A[i][j]*x[j]);
			}
			v[i] = x[i];
			x[i] = b[i] - soma;
		}
		NormaNum = NormaDen = 0;
		// Tira a norma infinita entre a iteração atual e a anterior
		for (i = 0; i < n; i++)
		{
			t = abs(x[i] - v[i]);
			if (t > NormaNum)
				NormaNum = t;
			if (abs(x[i]) > NormaDen)
				NormaDen = abs(x[i]);
		}
		NormaRel = NormaNum / NormaDen;
		imprimeIteracao(Iter,x,NormaRel);
		if (NormaRel <= Toler)
			break;
	}
	if (NormaRel <= Toler)
	{
		cout << "[+]Metodo encerrado corretamente" << endl;
		return (x);
	}
	else
	{
		cout << "[-]ERRO! Valor de tolerancia nao foi atingido!" << endl;
		return (NULL);
	}
}

double* calcX (double *x_ant, double lambda, double *d_ant)
{
	int i;
	// Variável auxiliar para preservar o valor de 'd_ant' após a operação com 'lambda'.
	double *aux = new double[n];
	for (i = 0; i < n; i++)
		aux[i] = lambda*d_ant[i];
	return (somaVetorVetor(x_ant,aux,n));
}

double* calcR (double **A, double *r_ant, double *d, double lambda)
{
	double *aux = multiplicaMatrizVetor(A,d,n);
	multiplicaEscalarVetor(aux,-lambda,n);
	return (somaVetorVetor(r_ant,aux,n));
}

double* calcD (double *r, double beta, double *d_ant)
{
	int i;	
	double *aux = new double[n];
	for (i = 0; i < n; i++)
		aux[i] = beta*d_ant[i];
	return (somaVetorVetor(r,aux,n));
	
}

double calcBeta (double *r, double *r_ant)
{
	double Num, Den;
	Num = multiplicaVetorVetor(r,r,n);
	Den = multiplicaVetorVetor(r_ant,r_ant,n);
	return (Num/Den);
}

double calcLambda (double **A, double *r, double *d)
{
	double Num, Den;	
	double *aux = multiplicaMatrizVetor(A,d,n);
	Num = multiplicaVetorVetor(r,r,n);
	Den = multiplicaVetorVetor(d,aux,n);
	return (Num/Den);
}

double* GradientesConjugados (double **A, double *b)
{
	double *d;					// Vetor de direções
	double lambda;					// Correção que tende a minimizar a função
	double *r;					// Vetor do resíduo, que deve diminuir a cada iteração
	double *r_anterior;				// Vetor do resíduo anterior
	double *x;					// Vetor solução
	double beta;					// Correção para pegar a próxima direção	
	double *aux;	
	int i;
	int Iter = 0;
	// Checar critério de convergência	
	if (!ehPositivaDefinida(A))
	{
		cout << "[-]ERRO! A matriz [A] nao eh positiva definida!" << endl;
		exit(1);
	}
	// Montando o vetor da solução inicial	
	x = new double[n];
	r_anterior = new double[n];
	for (i = 0; i < n; i++)
		x[i] = 0;
	aux = multiplicaMatrizVetor(A,x,n);
	multiplicaEscalarVetor(aux,-1,n);	
	r = somaVetorVetor(b,aux,n);
	d = somaVetorVetor(b,aux,n);

	while (Iter <= 10)
	{			
		lambda = calcLambda(A,r,d);
		x = calcX(x,lambda,d);
		copiaVetor(r,r_anterior,n);				// Fazer cópia do resíduo anterior
		r = calcR(A,r,d,lambda);			
		beta = calcBeta(r,r_anterior);
		d = calcD(r,beta,d);

		Iter++;
	}
	delete [] r;
	delete [] d;
	delete [] aux;

	cout << "[+]Algoritmo terminado com sucesso." << endl;
	imprimeSolucao(x);	
	return (x);	
}

// Lê a matriz dos coeficientes [A] do arquivo de entrada.
double** leMatriz ()
{
	int i, j;
	double **m;	
	fscanf(arq,"%d",&n);
	m = new double*[n];
	for (i = 0; i < n; i++)
	{
		m[i] = new double[n];
	}	
	
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			fscanf(arq,"%lf",&m[i][j]);
		}
	}
	return (m);
}

// Lê o vetor dos termos independentes {b} do arquivo de entrada.
double* leVetor ()
{	
	int i;	
	double *v = new double[n];
	for (i = 0; i < n; i++)
	{
		fscanf(arq,"%lf",&v[i]);
	}
	return (v);  
}

void imprimeSolucao (double *x)
{
	int i;
	cout << "\n=== SOLUCAO ===" << endl;
	for (i = 0; i < n; i++)
		cout << x[i] << endl;
}
