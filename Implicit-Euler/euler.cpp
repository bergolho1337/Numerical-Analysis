#include <math.h>
#include <iostream>
#define N 2
using namespace std;

void LU(double A[N][N],double B[N],double X[N]);

void newtonRaphson(double y[N],double yOld[N], double h, double t){
	double A[N][N];
	double F[N]; //-F
	double dY[N];
	double scalardY;


	cout<<"Vetor Y\n";
	for(int i=0;i<N;i++){
		cout<<y[i]<<"\n";
	}

	cout<<"Vetor YOld\n";
	for(int i=0;i<N;i++){
		cout<<yOld[i]<<"\n";
	}

	do{
		A[0][0]=1 - 2*h*y[0] - 2*h*y[1];
		A[0][1]=-2*h*y[0];
		A[1][0]=-h*t - 2*h*y[0]*sin(y[1]);
		A[1][1]=1 - h*y[0]*y[0]*cos(y[1]);

		F[0]=-y[0] + yOld[0] + h*y[0]*y[0] + 2*h*y[0]*y[1];
		F[1]=-y[1] + yOld[1] + h*t*y[0] + h*y[0]*y[0]*sin(y[1]);

		LU(A,F,dY);

		y[0]+=dY[0];
		y[1]+=dY[1];

		scalardY=sqrt(F[0]*F[0]+F[1]*F[1]);

		cout<<"Diferença escalar: "<<scalardY<<"\n";
		
	}while(scalardY>10e-10);

	cout<<"Vetor Y\n";
	for(int i=0;i<N;i++){
		cout<<y[i]<<"\n";
	}
}


void partialPivot(int *swap,double vet[N][N], int column);

void LU(double A[N][N],double B[N],double X[N]){

	int swap[N];
	for(int i=0;i<N;i++)swap[i]=i;

	double A_[N][N];

	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++)A_[i][j]=A[i][j];
	}
/*{
	cout<<"\nFuncão F\n";

	for(int i=0;i<N;i++){
		cout<<B[i]<<"\n";
	}
}*/
	for(int i=0;i<N;i++){
		partialPivot(swap,A,i);
		for(int j=i+1;j<N;j++){
			A[swap[j]][i]=A[swap[j]][i]/A[swap[i]][i];
			for(int k=i+1;k<N;k++){
				A[swap[j]][k]-=A[swap[j]][i]*A[swap[i]][k];
			}
		}
	}

	
	for(int i=0;i<N;i++)X[i]=B[swap[i]];

	for(int i=0;i<N;i++){
		for(int j=0;j<i+1;j++){
			if(i==j) X[i]=X[i];
			else X[i]-=A[swap[i]][j]*X[j];
		}
	}
	
	for(int i=N-1;i>=0;i--){
		for(int j=N-1;j>=i;j--){
			if(i==j) X[i]=X[i]/A[swap[i]][i];
			else X[i]-=A[swap[i]][j]*X[j];
		}
	}
/*{
	cout<<"Vetor Deslocamento\n";
	for(int i=0;i<N;i++){
		cout<<X[i]<<"\n";
	}

	double B_[N];
	for(int i=0;i<N;i++)B_[i]=0;

	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			B_[i]+=A_[i][j]*X[j];
		}
	}

	double error=0;
	for(int i=0;i<N;i++){
		error+=(B_[i]-B[i])*(B_[i]-B[i]);
	}

	error=sqrt(error);
	cout<<"\nErro: "<<scientific<<error<<"\n";
}*/
}



void partialPivot(int *swap,double vet[N][N], int column){
	double max=vet[column][column];
	if(max<0)max*=-1;
	int trade=column;
	for(int i=column+1;i<N;i++){
		double mod=vet[i][column];
		if(mod<0)mod*=-1;
		if(max<vet[i][column]){
			trade=i;
			max=mod;
		}
	}
	int aux=swap[column];
	swap[column]=swap[trade];
	swap[trade]=aux;
}



int main(){
	double y[2]={1,-1};
	double yOld[2]={1,-1};
	double t_0=0;


	double t=10;
	double h=1;


	for(double t_k=t_0+h ; t_k<=t ; t_k+=h){

		newtonRaphson(y,yOld,h,t_k);

		yOld[0]=y[0];
		yOld[1]=y[1];
	}

	cout<<"Vetor Y FINAL: \n";
	for(int i=0;i<N;i++){
		cout<<y[i]<<"\n";
	}

	return 0;
}
