#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>

using namespace std;

int nelem;
int n, m;

double* readRHS (const char rhs_filename[])
{
    ifstream in(rhs_filename);
    int i, j;
    double value;
    string str;

    getline(in,str);
    in >> n >> m;
    double *b = new double[n]();
    for (int i = 0; i < n; i++)
        in >> b[i];
    
    in.close();

    return b;
}

double* readMatrix (const char matrix_filename[])
{
    ifstream in(matrix_filename);
    int i, j;
    double value;
    string str;

    getline(in,str);
    in >> n >> m >> nelem;
    double *A = new double[nelem]();
    while (in >> i >> j >> value)
    {
        i--;
        j--;
        A[i*m+j] = value;
    }
    
    in.close();

    return A;
}

void printMatrix (const double *A)
{
    printf("\n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            printf("%.5lf ",A[i*n+j]);
        printf("\n");
    }
}

void printRHS (const double *b)
{
    printf("\n");
    for (int i = 0; i < n; i++)
    {
        printf("%.5lf\n",b[i]);
    }
}

void writeMatrixToLatex (const double *A)
{
    ofstream out("matrix_latex.txt");
    
    out << "\\begin{bmatrix}" << endl;
    for (int i = 0; i < n-1; i++)
    {
        for (int j = 0; j < n-1; j++)
            out << fixed << setprecision(5) << A[i*n+j] << " & ";
        out << fixed << setprecision(5) << A[i*n+(n-1)] << " \\\\" << endl;
    }
    for (int j = 0; j < n-1; j++)
            out << fixed << setprecision(5) << A[(n-1)*n+j] << " & ";
        out << fixed << setprecision(5) << A[(n-1)*n+(n-1)] << endl;
    out << "\\end{bmatrix}" << endl;
    
    out.close();
}

void writeRHSToLatex (const double *b)
{
    ofstream out("rhs_latex.txt");
    
    out << "\\begin{bmatrix}" << endl;
    for (int i = 0; i < n-1; i++)
    {
        out << fixed << setprecision(5) << b[i] << " \\\\" << endl;
    }
    out << fixed << setprecision(5) << b[n-1] << endl;
    out << "\\end{bmatrix}" << endl;
    
    out.close();
}

int main (int argc, char *argv[])
{
    double *A = readMatrix(argv[1]);
    double *b = readRHS(argv[2]);

    printMatrix(A);
    printRHS(b);

    writeMatrixToLatex(A);
    writeRHSToLatex(b);

    delete [] A;
    delete [] b;

    return 0;

}