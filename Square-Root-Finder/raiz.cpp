#include "raiz.h"

double avaliaFuncao (double x)
{
	//return ( ((pow(x+1,2))*(exp(x*x-2)))-1 );
	return (pow(x,3) + 3*x - 1);
}

void Bissecao (double a, double b)
{
	double x;
	int IterMax = 50, Iter;
	for (Iter = 0; Iter < IterMax; Iter++)
	{
		x = (a+b)/2;
		if (avaliaFuncao(a)*avaliaFuncao(x) < 0)
			b = x;
		else if (avaliaFuncao(a)*avaliaFuncao(x) > 0)
			a = x;
		else
		{
			cout << "[+] Raiz encontrada = " << x << endl;
			exit(1);
		}
	}
	cout << "[-] Numero maximo de iteracoes atingido." << endl;
	cout << "Aproximacao da raiz = " << x << endl;
	cout << "f(x) = " << avaliaFuncao(x) << endl;
}
