#include "linear-system-solver.h"

void swap (double **a, double **b)
{
    double *tmp = *a;
    *a = *b;
    *b = tmp;
}

double calcResidue (const double *A, const double *b, const double *x, const int n)
{
    double residue, sum;
    residue = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum = 0.0;
        for (int j = 0; j < n; j++)
            sum += A[i*n+j]*x[j];
        residue += pow(b[i] - sum,2);
    }
    return sqrt(residue);
}

int hasConverged (const double *A, const double *b, const double *x, const int n, const int iter)
{
    if (iter < MAX_ITER)
    {
        double residue = calcResidue(A,b,x,n);
        if (isnan(residue))
        {
            fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            fprintf(stderr,"[Linear-System-Solver] DANGER! Iter = %d || Residue = %.10lf\n",iter,residue);
            fprintf(stderr,"There is a problem with the method!\n");
            fprintf(stderr,"Check if the matrix is not singular!\n");
            fprintf(stderr,"If you are using Jacobi or Gauss_Seidel check if the matrix is diagonal dominant!\n");
            fprintf(stderr,"If you are using Conjugate Gradient check if the matrix is symmetric and positive definite!\n");
            fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            exit(EXIT_FAILURE);
        }
        if (residue < EPSILON)
            return 1;
        else
            return 0;
    }
    else
    {
        return 1; 
    }
}

void checkSolution (const double *A, const double *b, const double *x, const int n, const int iter)
{
    double residue = calcResidue(A,b,x,n);
    if (residue >= EPSILON)
    {
        fprintf(stderr,"*******************************************************************************************\n");
        fprintf(stderr,"[Linear-System-Solver] ERROR ! The method has not converged ! Iterations = %d\n",iter);
        fprintf(stderr,"*******************************************************************************************\n");
        exit(EXIT_FAILURE);
    }
    //fprintf(stdout,"[Linear-System-Solver] Number of iterations = %d\n",iter);
    //fprintf(stdout,"[Linear-System-Solver] Residue = %e\n",residue);
}

void Jacobi (double *A, double *b, const int n, double *x)
{
    int iter = 0;
    double sigma;
    double *x_aux = (double*)malloc(sizeof(double)*n);
    while (!hasConverged(A,b,x,n,iter))
    {
        // Find the next solution
        for (int i = 0; i < n; i++)
        {
            // Compute sigma
            sigma = 0.0;
            for (int j = 0; j < n; j++)
            {
                if (j != i)
                    sigma += A[i*n+j]*x[j];
            }
            // Calculate the new solution line
            x_aux[i] = (1/A[i*n+i])*(b[i] - sigma);
        }
        swap(&x,&x_aux);
        iter++;
    }
    checkSolution(A,b,x,n,iter);

    free(x_aux);
}

void Gauss_Seidel (double *A, double *b, const int n, double *x)
{
    int iter = 0;
    double sigma;
    while (!hasConverged(A,b,x,n,iter))
    {
        // Find the next solution
        for (int i = 0; i < n; i++)
        {
            // Compute sigma
            sigma = 0.0;
            for (int j = 0; j < n; j++)
            {
                if (j != i)
                    sigma += A[i*n+j]*x[j];
            }
            // Calculate the new solution line
            x[i] = (1/A[i*n+i])*(b[i] - sigma);
        }
        iter++;
    }
    checkSolution(A,b,x,n,iter);
}

// Conjugate gradient with/without a Jacobi preconditioner
void CG (double *A, double *b, const int n, double *x)
{
    int iter = 0;
    int use_jacobi;
    double rTr, rTz, pTAp, r1Tr1, r1Tz1;
    double alpha, beta;
    double error;
    double *r = (double*)malloc(sizeof(double)*n);
    double *z = (double*)malloc(sizeof(double)*n);
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);

    if (PRECONDITIONER) 
    {
        use_jacobi = 1;
        fprintf(stdout,"[Linear-System-Solver] Using Jacobi preconditioner\n");
    }
    else 
    {
        use_jacobi = 0;
    }               
    
    // Calculate r
    for (int i = 0; i < n; i++)
    {
        r[i] = 0.0;
        for (int j = 0; j < n; j++)
            r[i] += A[i*n+j]*x[j];
        r[i] = b[i] - r[i];
        p[i] = r[i];
    }

    while (!hasConverged(A,b,x,n,iter))
    {
        // Calc q
        pTAp = 0.0;
        rTr = 0.0;
        for (int i = 0; i < n; i++)
        {
            q[i] = 0.0;
            for (int j = 0; j < n; j++)
                q[i] += A[i*n+j]*p[j];
            pTAp += p[i] * q[i];
            rTr += r[i] * r[i];
        }
        
        // Calc alpha
        alpha = rTr / pTAp;


        // Update x and r
        for (int i = 0; i < n; i++)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
        }

        // Calc new r
        r1Tr1 = 0.0;
        for (int i = 0; i < n; i++)
            r1Tr1 += r[i] * r[i];

        // Calc beta
        beta = r1Tr1 / rTr;

        // Update p
        for (int i = 0; i < n; i++)
            p[i] = r[i] + beta * p[i];
        
        iter++;
    }
    checkSolution(A,b,x,n,iter);

    free(r);
    free(z);
    free(p);
    free(q);
}

// Biconjugate Gradient Method
void BiCG (double *A, double *b, const int n, double *x)
{
    int iter = 0;
    int use_jacobi;
    double rTr, rTz, pTAp, r1Tr1, r1Tz1;
    double alpha, beta;
    double error;
    double *x2 = (double*)malloc(sizeof(double)*n);
    double *r1 = (double*)malloc(sizeof(double)*n);
    double *r2 = (double*)malloc(sizeof(double)*n);
    double *p1 = (double*)malloc(sizeof(double)*n);
    double *p2 = (double*)malloc(sizeof(double)*n);
    double *q1 = (double*)malloc(sizeof(double)*n);
    double *q2 = (double*)malloc(sizeof(double)*n);
    
    // Set the second guess vector x2 (copy of the first one)
    for (int i = 0; i < n; i++)
        x2[i] = x[i];

    // Calculate r
    for (int i = 0; i < n; i++)
    {
        r1[i] = 0.0;
        r2[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            r1[i] += A[i*n+j]*x[j];
            r2[i] += x2[j]*A[j*n+i];
        }
            
        r1[i] = b[i] - r1[i];
        r2[i] = b[i] - r2[i];
        p1[i] = r1[i];
        p2[i] = r2[i];
    }

    while (!hasConverged(A,b,x,n,iter))
    {
        // Calc q1 and q2
        pTAp = 0.0;
        rTr = 0.0;
        for (int i = 0; i < n; i++)
        {
            q1[i] = 0.0;
            q2[i] = 0.0;
            for (int j = 0; j < n; j++)
            {
                q1[i] += A[i*n+j]*p1[j];
                q2[i] += A[j*n+i]*p2[j];
            }
            pTAp += p2[i] * q1[i];
            rTr += r1[i] * r2[i];
        }
        
        // Calc alpha
        alpha = rTr / pTAp;

        // Update x and r
        for (int i = 0; i < n; i++)
        {
            x[i] += alpha * p1[i];
            x2[i] += alpha * p2[i];
            r1[i] -= alpha * q1[i];
            r2[i] -= alpha * q2[i];
        }

        // Calc new residue
        r1Tr1 = 0.0;
        for (int i = 0; i < n; i++)
            r1Tr1 += r1[i] * r2[i];

        // Calc beta
        beta = r1Tr1 / rTr;

        // Update p
        for (int i = 0; i < n; i++)
        {
            p1[i] = r1[i] + beta * p1[i];
            p2[i] = r2[i] + beta * p2[i];
        }
        
        iter++;
    }
    checkSolution(A,b,x,n,iter);

    free(x2);
    free(r1);
    free(r2);
    free(p1);
    free(p2);
    free(q1);
    free(q2);
}