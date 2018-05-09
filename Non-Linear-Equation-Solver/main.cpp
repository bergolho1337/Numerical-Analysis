#include <iostream>
#include <cmath>

using namespace std;

// COLOQUE A FUNÇÃO QUE VC DESEJA OBTER A RAIZ AQUI !!!!
double avaliaFuncao (double x)
{
	//return ( 1 + (0.4 - 1)*x - 0.01*x*x*0.4 );
	return (-0.9472375 + 15.5925e-13*pow(x,4) - 0.9969145*x);
}

double avaliaDerivada (double x, double h)
{
	return ( (avaliaFuncao(x+h) - avaliaFuncao(x-h)) / (2*h) );  // Retorna um valor pequeno diferente de zero (pertubação).
}

double Newton (double x0, double h)
{
	double x_old, x_new, y_new;
	x_old = x0;
	x_new = x_old + 10*h;
	while (abs(x_new - x_old) > h)
	{
		y_new = avaliaFuncao(x_new);
		cout << "f(" << x_new << ") = " << y_new << endl;
		x_old = x_new;
		x_new = x_old - (y_new / avaliaDerivada(x_old,h)); 
	}
	return (x_new);
}

int main ()
{
	double s;	
	double x0 = 1;
	double h = 0.001;
	cout << "====== Resultado de Newton ========" << endl;
	s = Newton(x0,h);
	cout << "\nSolucao:\nx = " << s << endl;
}
