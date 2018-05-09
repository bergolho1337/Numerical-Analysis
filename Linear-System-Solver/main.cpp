#include <iostream>
#include "sistema.h"

// entrada.txt -> resolve-se pelo método LU (2.15)
// entrada2.txt -> resolve-se por Cholesky (a matriz é definida positiva) (2.23)
// entrada3.txt -> resolver por Gauss-Seidel (2.31)
// entrada4.txt -> resolver pelo método dos gradientes conjugados

using namespace std;

int main ()
{
	double *x;	
	double **A = leMatriz();
	double *b = leVetor();
	x = GradientesConjugados(A,b);
	//x = LU(A,b);
	//x = Cholesky(A,b);
	//x = Gauss_Seidel(A,b);
}

/* MODO DE USAR: no terminal */
// Compilar: $ make
// Rodar:    $ ./main
