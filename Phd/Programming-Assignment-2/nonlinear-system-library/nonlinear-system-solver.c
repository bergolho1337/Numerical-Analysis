#include "nonlinear-system-solver.h"

void NewtonFiniteDifference (double *J, double *f, const int n, double *x,\
                            set_linear_system_fn *solver_linear_system,\
                            set_problem_fn **functions)
{
    // TO DO: Implementar o metodo
    printTypeFiniteDifference();

    int iter = 0;
    // Allocate memory
    J = (double*)malloc(sizeof(double)*n*n);
    f = (double*)malloc(sizeof(double)*n);
    x = (double*)malloc(sizeof(double)*n);

    buildInitialGuess(x,n);
    do
    {
        if (iter % REBUILD_JACOBIAN == 0)
        {
            buildJacobian_FiniteDifferences(J,x,n,functions);
            /*
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    printf("%.10lf ",J[i*n+j]);
                printf("\n");
            }
            */
        }
            


        iter++;
    }while(!hasConverged(x,functions,n,iter));

    /*
    // DEBUG
    solver_linear_system(J,f,n,x);

    for (int i = 0; i < n; i++)
        printf("%lf\n",x[i]);
    
    printf("Testing problem functions\n");
    double v[3] = {0,0,0};
    for (int i = 0; i < 3; i++)
        printf("Function %d = %lf\n",i+1,functions[i](v));
    */

}

void buildJacobian_FiniteDifferences (double *J, double *x, const int n, set_problem_fn **functions)
{
    fprintf(stdout,"Rebuilding Jacobian matrix\n");
    double x_aux[n], x_aux2[n];
    double H12 = 0.5*H;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            memcpy(x_aux,x,sizeof(double)*n);
            switch (FINITE_DIFFERENCE)
            {
                // Forward Difference
                case 0: {
                            x_aux[j] = x[j] + H;

                            J[i*n+j] = (functions[i](x_aux) - functions[i](x)) / H;
                            break;
                        }
                // Central Difference
                case 1: {
                            memcpy(x_aux2,x,sizeof(double)*n);
                            x_aux[j] = x[j] + H12;
                            x_aux2[j] = x[j] - H12;

                            J[i*n+j] = (functions[i](x_aux) - functions[i](x_aux2)) / H; 
                            break;
                        }
            }
        }
    }

}

int hasConverged (const double *x, set_problem_fn **functions, const int n, const int iter)
{
    if (iter < MAX_ITER)
    {
        double residue = calcResidue(x,n,functions);
        if (isnan(residue))
        {
            fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            fprintf(stderr,"[Nonlinear-System-Solver] DANGER! Iter = %d\n",iter);
            fprintf(stderr,"There is a problem with the method!\n");
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

double calcResidue (const double *x, const int n, set_problem_fn **functions)
{
    double residue = 0.0;
    for (int i = 0; i < n; i++)
        residue += functions[i](x);
    return residue;
}

void buildInitialGuess (double *x, const int n)
{
    switch (GUESS_VECTOR)
    {
        // Null vector
        case 0: {
                    fprintf(stdout,"[Nonlinear-System-Solver] Building the null vector as initial guess\n");
                    for (int i = 0; i < n; i++)
                        x[i] = 0.0;
                    break;
                }
        // Random vector
        case 1: {
                    fprintf(stdout,"[Nonlinear-System-Solver] Building a random vector as initial guess\n");
                    srand(time(NULL));
                    for (int i = 0; i < n; i++)
                    {
                        if (rand() % 2 == 0)
                            x[i] = ((double)rand() / (double)RAND_MAX);
                        else
                            x[i] = -((double)rand() / (double)RAND_MAX);
                    }
                    break;
                }
        // Vector of one's where the even indexes are positive and the odds are negative
        case 2: {
                    fprintf(stdout,"[Nonlinear-System-Solver] Building a vector of one's as initial guess\n");
                    for (int i = 0; i < n; i++)
                    {
                        if (i % 2 == 0)
                            x[i] = 1.0;
                        else
                            x[i] = -1.0;
                    }
                    break;
                }
        // Custom vector
        case 3: {
                    fprintf(stdout,"[Nonlinear-System-Solver] Building a custom vector as initial guess\n");
                    fprintf(stdout,"Please enter the initial guess vector (n = %d):\n",n);
                    for (int i = 0; i < n; i++)
                    {
                        fprintf(stdout,"x[%d] = ",i);
                        scanf("%lf",&x[i]);
                    }
                    break;
                }
    }
}

void printTypeFiniteDifference ()
{
    switch (FINITE_DIFFERENCE)
    {
        case 0: fprintf(stdout,"[Nonlinear-System-Solver] Using Forward Finite Difference to approximate the derivatives\n");
                break;
        case 1: fprintf(stdout,"[Nonlinear-System-Solver] Using Central Finite Difference to approximate the derivatives\n");
                break;
    }
}