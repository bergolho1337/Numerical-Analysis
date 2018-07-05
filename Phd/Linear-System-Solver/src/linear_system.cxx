#include "../include/linear_system.h"

using namespace std;

void Usage (const char pname[])
{
    cout << PRINT_LINE << endl;
    cout << "Usage:> " << pname << " <matrix_file> <rhs_file> <id_method>" << endl;
    cout << "<matrix_file> = The matrix file should be on the MatrixMarket format" << endl;
    cout << "<rhs_file> = The RHS file should be on the MatrixMarket format" << endl;
    cout << "<id_method> = Identifier of the linear system solver" << endl;
    cout << "====================================================================================" << endl;
    cout << "**************************** DIRECT METHODS ****************************" << endl;
    cout << "\t1 - Forward Substitutions (Triangular)" << endl;
    cout << "\t2 - Backward Substitutions (Triangular)" << endl;
    cout << "\t3 - Gaussian Elimination" << endl;
    cout << "\t4 - LU Decomposition" << endl;
    cout << "************************** ITERACTIVE METHODS **************************" << endl;
    cout << "\t5 - Jacobi" << endl;
    cout << "\t6 - Gauss-Seidel" << endl;
    cout << "\t7 - Conjugate Gradient" << endl;
    cout << "\t8 - Biconjugate Gradient" << endl;
    cout << "\t9 - GMRES" << endl;
    cout << PRINT_LINE << endl;
    cout << "Example: " << pname << " Inputs/Square/Non-Symmetric/e05r0000.mtx " << \
                                      "Inputs/Square/Non-Symmetric/e05r0000_rhs1.mtx 8" << endl;
    cout << PRINT_LINE << endl;
    
}

LinearSystem::LinearSystem (int argc, char *argv[])
{
    // Read matrix
    readMatrix(argv[1]);
    //printMatrix();
    
    // Read RHS
    readRHS(argv[2]);
    //printRHS();

    // Read type of the linear system solver
    method = atoi(argv[3]);
}

LinearSystem::~LinearSystem ()
{
    cout << "[+] Freeing memory ..." << endl;
    free(A[0]);
    free(A);
    free(b);
    free(x);
    // Remembering to free the LU matrix ...
    if (method == 4) 
    {
        free(LU[0]);
        free(LU);
    }
}

void LinearSystem::readMatrix (const char fname[])
{
    string str;
    int n, m, k;
    int lin, col;
    double val;

    ifstream in(fname);
    getline(in,str);
    in >> n >> m >> k;
    N = n;
    A = allocMatrix(N,N);
    for (int i = 0; i < k; i++)
    {
        in >> lin >> col >> val;
        A[lin-1][col-1] = val;
    }
    in.close();
}

void LinearSystem::printMatrix ()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
}

// The addresses of the lines are contigous on memory
double** LinearSystem::allocMatrix (const int n, const int m)
{
    cout << "[+] Allocating memory for matrix ..." << endl;
    double **A = (double**)malloc(sizeof(double*)*n);
    A[0] = (double*)malloc(sizeof(double)*n*m);
    // Hacking the addresses of the lines
    for (int i = 1; i < n; i++)
        A[i] = A[0] + m*i;
    return A;
}

void LinearSystem::readRHS (const char fname[])
{
    string str;
    int n, k;
    double val;

    ifstream in(fname);
    getline(in,str);
    in >> n >> k;
    allocRHS(N);
    for (int i = 0; i < N; i++)
    {
        in >> val;
        b[i] = val;
    }
    in.close();
}

void LinearSystem::allocRHS (const int n)
{
    cout << "[+] Allocating memory for right-hand side ..." << endl;
    b = (double*)malloc(sizeof(double)*n);
}

double* LinearSystem::allocVector (const int n)
{
    cout << "[+] Allocating memory for vector ..." << endl;
    double *v = (double*)malloc(sizeof(double)*n);
    return v;
}

void LinearSystem::printRHS ()
{
    for (int i = 0; i < N; i++)
        cout << b[i] << endl;
}

void LinearSystem::solve ()
{
    // Initial guess is a zero vector
    x = (double*)calloc(N,sizeof(double));
    // Initial guess is the 'rhs' vector
    for (int i = 0; i < N; i++)
        x[i] = b[i];

    switch (method)
    {
        case 1: ForwardSubstitution();
                break;
        case 2: BackwardSubstitution();
                break;
        case 3: GaussianElimination();
                break;
        case 4: LUDecomposition();
                break;
        case 5: Jacobi();
                break;
        case 6: Gauss_Seidel();
                break;
        case 7: CG();
                break;
        case 8: BiCG();
                break;
        case 9: GMRES();
                break;
    }

    checkSolution();
}

// Solves an Inferior Triangular linear system
void LinearSystem::ForwardSubstitution ()
{
    cout << "[+] Solving linear system using Forward Substitution ..." << endl;
    x[0] = b[0] / A[0][0];
    for (int i = 1; i < N; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < i; j++)
            sum += A[i][j] * x[j];
        x[i] = (b[i] - sum) / A[i][i];
    }
}

// Solves an Superior Triangular linear system
void LinearSystem::BackwardSubstitution ()
{
    cout << "[+] Solving linear system using Backward Substitution ..." << endl;
    x[N-1] = b[N-1] / A[N-1][N-1];
    for (int i = N-2; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i+1; j < N; j++)
            sum += A[i][j] * x[j];
        x[i] = (b[i] - sum) / A[i][i];
    }
}

void LinearSystem::reduceToRowEchelonForm ()
{
    // For each line
    for (int i = 0; i < N-1; i++)
    {
        for (int j = i+1; j < N; j++)
        {
            double m = A[j][i] / A[i][i];
            // Multiply all the line using the multiplier 'm'
            for (int k = i; k < N; k++)
            {
                A[j][k] -= (m * A[i][k]);
            }
            b[j] -= m * b[i];
        }
    }
}

// Solves a general linear system without partial pivoting
void LinearSystem::GaussianElimination ()
{
    cout << "[+] Solving linear system using Gaussian Elimination ..." << endl;
    // A[0][0] is the pivot element
    // TO DO: include partial pivoting

    reduceToRowEchelonForm();

    // Solve the resulting system using backward substitution
    x[N-1] = b[N-1] / A[N-1][N-1];
    for (int i = N-2; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i+1; j < N; j++)
            sum += A[i][j] * x[j];
        x[i] = (b[i] - sum) / A[i][i];
    }
}

void LinearSystem::copyMatrixLU ()
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            LU[i][j] = A[i][j];
}

void LinearSystem::choosePivot (int &pivot_line, double &Amax, const int i)
{
    pivot_line = i;
    Amax = fabs(LU[i][i]);

    // For all the lines below 'i' search for the maximum value in module
    for (int j = i+1; j < N; j++)
    {
        double value = fabs(LU[j][i]);
        if (value > Amax)
        {
            Amax = value;
            pivot_line = j; 
        }
    }
}

void LinearSystem::switchLines (int pivot[], const int pivot_line, const int i)
{
    // Copy the role line
    // TO DO: maybe a pointer swap will be faster ...
    for (int j = 0; j < N; j++)
    {
        double tmp = LU[i][j];
        LU[i][j] = LU[pivot_line][j];
        LU[pivot_line][j] = tmp;
    }
    int m = pivot[i];
    pivot[i] = pivot[pivot_line];
    pivot[pivot_line] = m;
}

void LinearSystem::LUDecomposition ()
{
    cout << "[+] Solving linear system using LU Decomposition ..." << endl;

    // Copy the matrix A to LU
    LU = allocMatrix(N,N);
    copyMatrixLU();
    
    // Auxiliary vector y and pivot, they can be static arrays ...
    double y[N];
    int pivot[N];

    int pivot_line;
    double Amax;

    for (int i = 0; i < N; i++)
        pivot[i] = i;

    // For each line
    for (int i = 0; i < N-1; i++)
    {
        choosePivot(pivot_line,Amax,i);
        // The pivot element changed ? If so we need to switch lines
        if (pivot_line != i)
        {
            switchLines(pivot,pivot_line,i);
        }

        // Check if the matrix is not singular
        // If not, apply a Gaussian Elimination on the current line
        if (fabs(LU[i][i]) != 0.0)
        {
            double r = 1 / LU[i][i];
            for (int j = i+1; j < N; j++)
            {
                double m = LU[j][i] * r;

                // Write the multiplier on the inferior triangular matrix L
                LU[j][i] = m;

                // Write the results on the superior triangular matrix U
                for (int k = i+1; k < N; k++)
                {
                    LU[j][k] -= (m * LU[i][k]);
                }
            }
        }
    }

    // Forward substitution using the pivot array    
    int k = pivot[0];
    y[0] = b[k];
    for (int i = 1; i < N; i++)
    {
        double sum = 0.0;
        for (int j = 0; j <= i; j++)
            sum += LU[i][j] * y[j];
        k = pivot[i];
        y[i] = (b[k] - sum);
    }

    // Backward substitution
    x[N-1] = y[N-1] / LU[N-1][N-1];
    for (int i = N-2; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i+1; j < N; j++)
            sum += LU[i][j] * x[j];
        x[i] = (y[i] - sum) / LU[i][i];
    }

}

void LinearSystem::Jacobi ()
{
    cout << "[+] Solving linear system using Jacobi ..." << endl;
    int iter = 0;
    double sigma;
    double *x_aux = (double*)malloc(sizeof(double)*N);
    while (!hasConverged(iter))
    {
        // Find the next solution
        for (int i = 0; i < N; i++)
        {
            // Compute sigma
            sigma = 0.0;
            for (int j = 0; j < N; j++)
            {
                if (j != i)
                    sigma += A[i][j]*x[j];
            }
            // Calculate the new solution line
            x_aux[i] = (1/A[i][i])*(b[i] - sigma);
        }
        swap(x,x_aux);
        iter++;
    }
    free(x_aux);
    cout << "Number of iterations = " << iter << endl;
}

void LinearSystem::Gauss_Seidel ()
{
    cout << "[+] Solving linear system using Gauss-Seidel ..." << endl;
    int iter = 0;
    double sigma;
    while (!hasConverged(iter))
    {
        // Find the next solution
        for (int i = 0; i < N; i++)
        {
            // Compute sigma
            sigma = 0.0;
            for (int j = 0; j < N; j++)
            {
                if (j != i)
                    sigma += A[i][j]*x[j];
            }
            // Calculate the new solution line
            x[i] = (1/A[i][i])*(b[i] - sigma);
        }
        iter++;
    }
    cout << "Number of iterations = " << iter << endl;
}

void LinearSystem::CG ()
{
    cout << "[+] Solving linear system using Conjugate Gradient ..." << endl;
    const double TOL = 1.0e-16;
    int iter = 0;
    bool use_jacobi = false;
    double rTr, rTz, pTAp, r1Tr1, r1Tz1;
    double alpha, beta;
    double error;
    double *r = (double*)malloc(sizeof(double)*N);
    double *z = (double*)malloc(sizeof(double)*N);
    double *p = (double*)malloc(sizeof(double)*N);
    double *q = (double*)malloc(sizeof(double)*N);
    
    // Calculate r
    for (int i = 0; i < N; i++)
    {
        r[i] = 0.0;
        for (int j = 0; j < N; j++)
            r[i] += A[i][j]*x[j];
        r[i] = b[i] - r[i];
        p[i] = r[i];
    }

    while (!hasConverged(iter))
    {
        // Calc q
        pTAp = 0.0;
        rTr = 0.0;
        for (int i = 0; i < N; i++)
        {
            q[i] = 0.0;
            for (int j = 0; j < N; j++)
                q[i] += A[i][j]*p[j];
            pTAp += p[i] * q[i];
            rTr += r[i] * r[i];
        }
        
        // Calc alpha
        alpha = rTr / pTAp;


        // Update x and r
        for (int i = 0; i < N; i++)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
        }

        // Calc new r
        r1Tr1 = 0.0;
        for (int i = 0; i < N; i++)
            r1Tr1 += r[i] * r[i];

        // Calc beta
        beta = r1Tr1 / rTr;

        // Update p
        for (int i = 0; i < N; i++)
            p[i] = r[i] + beta * p[i];
        
        iter++;
    }
    cout << "Number of iterations = " << iter << endl;

    free(r);
    free(z);
    free(p);
    free(q);
}

void LinearSystem::BiCG ()
{
    cout << "[+] Solving linear system using Biconjugate Gradient ..." << endl;
    const double TOL = 1.0e-16;
    int iter = 0;
    bool use_jacobi = false;
    double rTr, rTz, pTAp, r1Tr1, r1Tz1;
    double alpha, beta;
    double error;
    double *x2 = (double*)malloc(sizeof(double)*N);
    double *r1 = (double*)malloc(sizeof(double)*N);
    double *r2 = (double*)malloc(sizeof(double)*N);
    double *p1 = (double*)malloc(sizeof(double)*N);
    double *p2 = (double*)malloc(sizeof(double)*N);
    double *q1 = (double*)malloc(sizeof(double)*N);
    double *q2 = (double*)malloc(sizeof(double)*N);
    //double *z = (double*)malloc(sizeof(double)*N);
    
    // Set the second guess vector x2 (copy of the first one)
    for (int i = 0; i < N; i++)
        x2[i] = x[i];

    // Calculate r
    for (int i = 0; i < N; i++)
    {
        r1[i] = 0.0;
        r2[i] = 0.0;
        for (int j = 0; j < N; j++)
        {
            r1[i] += A[i][j]*x[j];
            r2[i] += x2[j]*A[j][i];
        }
            
        r1[i] = b[i] - r1[i];
        r2[i] = b[i] - r2[i];
        p1[i] = r1[i];
        p2[i] = r2[i];
    }

    while (!hasConverged(iter))
    {
        // Calc q1 and q2
        pTAp = 0.0;
        rTr = 0.0;
        for (int i = 0; i < N; i++)
        {
            q1[i] = 0.0;
            q2[i] = 0.0;
            for (int j = 0; j < N; j++)
            {
                q1[i] += A[i][j]*p1[j];
                q2[i] += A[j][i]*p2[j];
            }
            pTAp += p2[i] * q1[i];
            rTr += r1[i] * r2[i];
        }
        
        // Calc alpha
        alpha = rTr / pTAp;

        // Update x and r
        for (int i = 0; i < N; i++)
        {
            x[i] += alpha * p1[i];
            x2[i] += alpha * p2[i];
            r1[i] -= alpha * q1[i];
            r2[i] -= alpha * q2[i];
        }

        // Calc new residue
        r1Tr1 = 0.0;
        for (int i = 0; i < N; i++)
            r1Tr1 += r1[i] * r2[i];

        // Calc beta
        beta = r1Tr1 / rTr;

        // Update p
        for (int i = 0; i < N; i++)
        {
            p1[i] = r1[i] + beta * p1[i];
            p2[i] = r2[i] + beta * p2[i];
        }
        
        iter++;
    }
    cout << "Number of iterations = " << iter << endl;

    free(x2);
    free(r1);
    free(r2);
    free(p1);
    free(p2);
    free(q1);
    free(q2);
}

void LinearSystem::GMRES ()
{
    cout << "[+] Solving linear system using GMRES ..." << endl;
}

void LinearSystem::checkSolution ()
{
    double residue = calcResidue();
    cout << "Norm of the residue = " << setprecision(20) << residue << endl;
    cout << PRINT_LINE << endl;
    cout << "Solution" << endl;
    for (int i = 0; i < N; i++)
        cout << setprecision(20) << x[i] << endl;
    cout << PRINT_LINE << endl;
}

bool LinearSystem::hasConverged (const int iter)
{
    if (iter < MAX_ITER)
    {
        double residue = calcResidue();
        if (isnan(residue))
        {
            cout << "[-] DANGER! Iter = " << iter << endl;
            cout << "Problem with the method!" << endl;
            exit(EXIT_FAILURE);
        }
        if (residue < EPSILON)
            return true;
        else
            return false;
    }
    else
    {
        return true; 
    }
}

double LinearSystem::calcResidue ()
{
    double residue, sum;
    residue = 0.0;
    for (int i = 0; i < N; i++)
    {
        sum = 0.0;
        for (int j = 0; j < N; j++)
            sum += A[i][j]*x[j];
        residue += pow(b[i] - sum,2);
    }
    return sqrt(residue);
}

void print (const char name[], double *v, const int N)
{
    cout << name << endl;
    for (int i = 0; i < N; i++)
        cout << v[i] << endl;
}