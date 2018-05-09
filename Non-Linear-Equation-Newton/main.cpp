#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>

using namespace std;

double avaliaFuncao1 (double t, double *y)
{
	return (-pow(y[0],2) - y[1]);
}

double avaliaFuncao2 (double t, double *y)
{
	return (2*y[0] - pow(y[1],3));
}

void partialPivot(int *swap,double **vet, int column, int N){

    int i;
	double max=vet[column][column];
	int trade=column;

	for(i=column+1;i<N;i++){

		double mod=vet[i][column];
		if(mod<0)mod*=-1;
		if(max<vet[i][column]){
			trade=i;
			max=vet[i][column];
		}
	}
	int aux=swap[column];
	swap[column]=swap[trade];
	swap[trade]=aux;
}

double* LU(double **A,double *B, int N){

    int i,j,k;
    int swap[N];

    for(i=0;i<N;i++)swap[i]=i;


	for(i=0;i<N;i++){
		partialPivot(swap,A,i,N);
		for(j=i+1;j<N;j++){
			A[swap[j]][i]=A[swap[j]][i]/A[swap[i]][i];
			for(k=i+1;k<N;k++){
				A[swap[j]][k]-=A[swap[j]][i]*A[swap[i]][k];
			}
		}
	}
/*
    printf("Matriz LU\n");
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			printf("%f \t",A[swap[i]][j]);
		}
		printf("\n");
	}
*/
    double Y[N];
	for(i=0;i<N;i++)Y[i]=B[swap[i]];
	double *X = (double*)malloc(N*sizeof(double));
	for(i=0;i<N;i++){
		for(j=0;j<i+1;j++){
			if(i==j) Y[i]=Y[i];
			else Y[i]-=A[swap[i]][j]*Y[j];
		}
	}

/*
	printf("\nVetor Y\n");
	for(i=0;i<N;i++){
		printf("%f",Y[i]);
		printf("\n");
	}

*/
	for(i=0;i<N;i++)X[i]=Y[i];
	for(i=N-1;i>=0;i--){
		for(j=N-1;j>=i;j--){
			if(i==j) X[i]=X[i]/A[swap[i]][i];
			else X[i]-=A[swap[i]][j]*X[j];
		}
	}
/*
	printf("\nVetor X\n");
	for(i=0;i<N;i++){
		printf("%f",X[i]);
		printf("\n");
	}
*/
    return X;


}


double* NewtonRaphson (double *y_old, double h, double t)
{
	int i;	
	double **J;
	double *F, *H, *y_new;
	double y_first[2];
	double dif[2];
	double y1, y2;
	double epsilon = 0.01;
	double Toler = 1.0e-05;
	y_new = new double[2];
	F = new double[2];
	J = new double*[2];
	for (i = 0; i < 2; i++)
		J[i] = new double[2];
		
	y_first[0] = y_old[0];
	y_first[1] = y_old[1];	
	do
	{
			
		// Monta vetor termos independentes
		F[0] = -1*(y_first[0] - y_old[0] + (h*avaliaFuncao1(t,y_old)));
		F[1] = -1*(y_first[1] - y_old[1] + (h*avaliaFuncao2(t,y_old)));
		// Monta matriz jacobiana
		dif[0] = y_old[0] + epsilon;
		dif[1] = y_old[1];
		J[0][0] = ((y_first[0] + dif[0] + (h*avaliaFuncao1(t,dif))) - (y_first[0] + y_old[0] + (h*avaliaFuncao1(t,y_old)))) / epsilon - 1;
		J[1][0] = ((y_first[1] + dif[1] + (h*avaliaFuncao2(t,dif))) - (y_first[1] + y_old[1] + (h*avaliaFuncao2(t,y_old)))) / epsilon;
		dif[0] = y_old[0] - epsilon;
		dif[1] = y_old[1] + epsilon;
		J[0][1] = ((y_first[0] + dif[0] + (h*avaliaFuncao1(t,dif))) - (y_first[0] + y_old[0] + (h*avaliaFuncao1(t,y_old)))) / epsilon;
		J[1][1] = ((y_first[1] + dif[1] + (h*avaliaFuncao2(t,dif))) - (y_first[1] + y_old[1] + (h*avaliaFuncao2(t,y_old)))) / epsilon - 1;

		H = LU(J,F,2);
		y_old[0] = y_old[0] + H[0];
		y_old[1] = y_old[1] + H[1];

		y1 = y_first[0] - y_old[0] + (h*avaliaFuncao1(t,y_old));
		y2 = y_first[1] - y_old[1] + (h*avaliaFuncao2(t,y_old));

		printf("%e\n", sqrt(pow(y1,2) + pow(y2,2)));
		
	} while ( sqrt(pow(y1,2) + pow(y2,2)) > Toler ); // Verifica se a solução está próxima do zero da função
		
	cout << "\nSolucao ----- t = " << t << endl;
	printf("%e", sqrt(pow(y1,2) + pow(y2,2)));
 	delete [] J;
	delete [] F;
	y_new[0] = y_old[0];
	y_new[1] = y_old[1];
	y_old[0] = y_first[0];
	y_old[1] = y_first[1];
	return (y_new);
}	

void geraGrafico ()
{
	ofstream arq;
	int tick;
	arq.open("graph.plt");
	arq << "set title \"Método Euler Implícito\"" << endl;
	arq << "set grid" << endl;
	arq << "set xlabel \"t\"" << endl;
	arq << "set ylabel \"y\"" << endl;
	arq << "set terminal png" << endl;
	arq << "set output \'sistema.png\'" << endl;
	arq << "plot \'dados.dat\' using 1:2 title \"y1\" w l, \'dados.dat\' using 1:3 title \"y2\" w l" << endl;
	arq.close();	
	
	tick = system("gnuplot graph.plt");
	if (!tick)
		cout << "OK" << endl;
	else
		cout << "[-] ERRO! Arquivo graph.plt com erro." << endl;
}

int main ()
{
	ofstream arq;	
	double y0[2] = {1,1};
	double f[2];
	double *y_new;
	double t0 = 0, tmax = 20, h = 0.1, t;
	arq.open("dados.dat");
	for (t = t0; t < tmax; t += h)
	{
		arq << t << "\t" << y0[0] << "\t" << y0[1] << endl;
		// Descobrir y-j+1
		y_new = NewtonRaphson(y0,h,t);
		f[0] = avaliaFuncao1(t+h,y_new);
		f[1] = avaliaFuncao2(t+h,y_new);
		y_new[0] = y0[0] + h*f[0];
		y_new[1] = y0[1] + h*f[1];
		y0[0] = y_new[0];
		y0[1] = y_new[1];
	}
	arq.close();
	geraGrafico();
}
