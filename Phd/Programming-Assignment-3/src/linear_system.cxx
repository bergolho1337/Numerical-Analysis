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
    cout << "\t3 - Gaussian Elimination no Pivoting" << endl;
    cout << "\t4 - LU Decomposition with Partial Pivoting" << endl;
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
    // Read type of the linear system solver
    method = atoi(argv[3]);

    // Read matrix
    readMatrix(argv[1]);
    //printMatrix();
    
    // Read RHS
    readRHS(argv[2]);
    //printRHS();
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
    if (method == 4) 
        LU = allocMatrix(N,N);
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
            cout << fixed << setprecision(10) << A[i][j] << " ";
        cout << endl;
    }
}

// The addresses of the lines are contigous on memory
double** allocMatrix (const int n, const int m)
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
    b = allocVector(N);
    for (int i = 0; i < N; i++)
    {
        in >> val;
        b[i] = val;
    }
    in.close();
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
        cout << fixed << setprecision(10) << b[i] << endl;
}

void LinearSystem::write ()
{
    ofstream out("scripts/error/solution.dat");
    for (int i = 0; i < N; i++)
        out << fixed << setprecision(10) << x[i] << endl;
    out.close();
}

void LinearSystem::cond ()
{
    cout << "[+] Calculating the inverse of A ..." << endl;
    double **Ainv = allocMatrix(N,N);
    double *e = allocVector(N);
    double *v = allocVector(N);
    // Calculate each column of the inverse matriz by solving a linear system
    for (int k = 0; k < N; k++)
    {
        cout << "[!] Calculting column " << k << " ..." << endl;
        
        buildElementaryVector(e,k,N);
        CG(A,e,v,N);
        insertColumnMatrix(Ainv,v,k,N);
        
        cout << endl;
    }

    // Test the inverse
    cout << "[+] Testing A inverse ..." << endl;
    MatrixMatrixMultiplication(A,Ainv,N);

    // Calculate the infinite matrix norm of A and Ainv
    double normA = calcMatrixNormInf(A,N);
    double normAinv = calcMatrixNormInf(Ainv,N);

    cout << PRINT_DANGER << endl;
    cout << fixed << setprecision(10) << "Condition number = " << normA*normAinv << endl;
    cout << PRINT_DANGER << endl;

    // Free allocated memory 
    free(Ainv[0]);
    free(Ainv);
    free(e);
    free(v);    
}

void LinearSystem::solve ()
{
    // Initial guess is a zero vector
    x = allocVector(N);
    // Initial guess is the 'rhs' vector
    for (int i = 0; i < N; i++)
        //x[i] = 101.0;
        x[i] = b[i];

    switch (method)
    {
        case 1: ForwardSubstitution(A,b,x,N);
                break;
        case 2: BackwardSubstitution(A,b,x,N);
                break;
        case 3: GaussianElimination(A,b,x,N);
                break;
        case 4: LUDecomposition(A,LU,b,x,N);
                break;
        case 5: Jacobi(A,b,x,N);
                break;
        case 6: Gauss_Seidel(A,b,x,N);
                break;
        case 7: CG(A,b,x,N);
                break;
        case 8: BiCG(A,b,x,N);
                break;
        case 9: GMRES(A,b,x,N);
                break;
    }

    if (!checkSolution(A,x,b,N))
        refineSolution(A,x,b,N);

}

// Solves an Inferior Triangular linear system
void ForwardSubstitution (double **A, double *b, double *x, const int N)
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
void BackwardSubstitution (double **A, double *b, double *x, const int N)
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

void reduceToRowEchelonForm (double **A, double *b, const int N)
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

// Solves a general linear system WITHOUT partial pivoting
void GaussianElimination (double **A, double *b, double *x, const int N)
{
    cout << "[+] Solving linear system using Gaussian Elimination ..." << endl;
    // A[0][0] is the pivot element
    // TO DO: include partial pivoting

    reduceToRowEchelonForm(A,b,N);

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

void copyMatrixLU (double **A, double **LU, const int N)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            LU[i][j] = A[i][j];
}

void choosePivot (double **LU, int &pivot_line, double &Amax, const int i, const int N)
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

void switchLines (double **LU, int pivot[], const int pivot_line, const int i, const int N)
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

// LU Decomposition with partial pivoting
void LUDecomposition (double **A, double **LU, double *b, double *x, const int N)
{
    cout << "[+] Solving linear system using LU Decomposition with partial pivoting ..." << endl;

    // Copy the matrix A to LU
    copyMatrixLU(A,LU,N);
    
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
        choosePivot(LU,pivot_line,Amax,i,N);
        // The pivot element changed ? If so we need to switch lines
        if (pivot_line != i)
        {
            switchLines(LU,pivot,pivot_line,i,N);
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

void Jacobi (double **A, double *b, double *x, const int N)
{
    cout << "[+] Solving linear system using Jacobi ..." << endl;
    int iter = 0;
    double sigma;
    double *x_aux = (double*)malloc(sizeof(double)*N);
    while (!hasConverged(A,x,b,iter,N))
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
    if (calcResidue(A,x,b,N) >= EPSILON)
    {
        cout << PRINT_DANGER << endl;
        cout << "[-] ERROR ! The method has not converged !" << endl;
        cout << PRINT_DANGER << endl;
    }
    cout << "Number of iterations = " << iter << endl;
}

void Gauss_Seidel (double **A, double *b, double *x, const int N)
{
    cout << "[+] Solving linear system using Gauss-Seidel ..." << endl;
    int iter = 0;
    double sigma;
    while (!hasConverged(A,x,b,iter,N))
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
    if (calcResidue(A,x,b,N) >= EPSILON)
    {
        cout << PRINT_DANGER << endl;
        cout << "[-] ERROR ! The method has not converged !" << endl;
        cout << PRINT_DANGER << endl;
    }
        
    cout << "Number of iterations = " << iter << endl;
}

void CG (double **A, double *b, double *x, const int N)
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

    while (!hasConverged(A,x,b,iter,N))
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

void BiCG (double **A, double *b, double *x, const int N)
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

    while (!hasConverged(A,x,b,iter,N))
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

void GMRES (double **A, double *b, double *x, const int N)
{
    cout << "[+] Solving linear system using GMRES ..." << endl;
}

void buildElementaryVector (double *e, const int k, const int N)
{
    for (int i = 0; i < N; i++)
    {
        if (i == k)
            e[i] = 1.0;
        else
            e[i] = 0.0;
    }
}

// Copy the vector 'v' into the 'k' column of 'A'
void insertColumnMatrix (double **A, double *v, const int k, const int N)
{
    for (int i = 0; i < N; i++)
        A[i][k] = v[i];
}

bool checkSolution (double **A, const double *x, const double *b, const int N)
{
    double residue = calcResidue(A,x,b,N);
    cout << "Norm of the residue = " << setprecision(20) << residue << endl;
    cout << PRINT_LINE << endl;
    cout << "Solution" << endl;
    for (int i = 0; i < N; i++)
        cout << setprecision(20) << x[i] << endl;
    cout << PRINT_LINE << endl;

    // If the residue is greater than the tolerance we will apply a refinement
    return (residue > EPSILON) ? false : true;
}

bool hasConverged (double **A, const double *x, const double *b, const int iter, const int N)
{
    if (iter < MAX_ITER)
    {
        double residue = calcResidue(A,x,b,N);
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

// Using Norm-2
double calcResidue (double **A, const double *x, const double *b, const int N)
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

// Using Norm-INF
double calcResidueInf (double **A, const double *x, const double *b, const int N)
{
    double norm = __DBL_MIN__;
    for (int i = 0; i < N; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < N; j++)
            sum += A[i][j]*x[j];
        double residue = fabs(b[i] - sum);
        if (residue > norm)
            norm = residue;
    }
    return norm;
}

void compResidue (double **A, const double *x, const double *b, double *r, const int N)
{
    double sum;
    for (int i = 0; i < N; i++)
    {
        sum = 0.0;
        for (int j = 0; j < N; j++)
        {
            sum += A[i][j]*x[j];
        }
        r[i] = b[i] - sum;
    }
}

bool refineSolution (double **A, double *x, double *b, const int N)
{
    cout << "[+] Refining solution ..." << endl;

    int iter = 0;
    double r[N], d[N];
    
    compResidue(A,x,b,r,N);
        
        Gauss_Seidel(A,r,d,N);
        
        for (int i = 0; i < N; i++)
            x[i] += d[i];

    bool ok = checkSolution(A,x,b,N);

    return ok;
}

void MatrixMatrixMultiplication (double **A, double **B, const int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double value = 0.0;
            for (int k = 0; k < N; k++)
            {
                value += A[i][k]*B[k][j];
            }
            cout << fixed << setprecision(10) << value << " "; 
        }
        cout << endl;
    }
}

// Calculate the infinite norm of a matrix A
double calcMatrixNormInf (double **A, const int N)
{
    double norm = 0.0;
    for (int i = 0; i < N; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < N; j++)
            sum += fabs(A[i][j]);
        norm = max(norm,sum);
    }
    return norm;
}

void print (const char name[], double *v, const int N)
{
    cout << name << endl;
    for (int i = 0; i < N; i++)
        cout << v[i] << endl;
}