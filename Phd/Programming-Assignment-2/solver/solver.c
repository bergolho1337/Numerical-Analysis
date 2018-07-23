#include "solver.h"

struct solver_data* new_solver_data (int argc, char *argv[])
{
    struct solver_data *s = (struct solver_data*)malloc(sizeof(struct solver_data));
    s->problem_id = atoi(argv[1]);
    s->linear_system_solver_id = atoi(argv[2]);
    s->nonlinear_system_solver_id = atoi(argv[3]);

    s->problem = new_problem(s->problem_id);
    //printProblem(s->problem);     // DEBUG

    s->linear_system_solver = new_linear_system(s->linear_system_solver_id);

    /*
    DEBUG
    s->linear_system_solver->n = 4;

    s->linear_system_solver->A = (double*)malloc(sizeof(double)*4*4);
    s->linear_system_solver->A[0] = 4.000000000000;
    s->linear_system_solver->A[1] = -1.000000000000;
    s->linear_system_solver->A[2] = 0.000000000000;
    s->linear_system_solver->A[3] = -1.000000000000;
    s->linear_system_solver->A[4] = 1.000000000000;
    s->linear_system_solver->A[5] = -2.000000000000;
    s->linear_system_solver->A[6] = 1.000000000000;
    s->linear_system_solver->A[7] = 0.000000000000;
    s->linear_system_solver->A[8] = 0.000000000000;
    s->linear_system_solver->A[9] = 4.000000000000;
    s->linear_system_solver->A[10] = -4.000000000000;
    s->linear_system_solver->A[11] = 1.000000000000;
    s->linear_system_solver->A[12] = 5.000000000000;
    s->linear_system_solver->A[13] = 0.000000000000;
    s->linear_system_solver->A[14] = 5.000000000000;
    s->linear_system_solver->A[15] = -1.000000000000;

    s->linear_system_solver->b = (double*)malloc(sizeof(double)*4);
    s->linear_system_solver->b[0] = 1.000000000000;
    s->linear_system_solver->b[1] = -2.000000000000;
    s->linear_system_solver->b[2] = -3.000000000000;
    s->linear_system_solver->b[3] = 4.000000000000;

    s->linear_system_solver->x = (double*)malloc(sizeof(double)*4);
    s->linear_system_solver->solver(s->linear_system_solver->A,s->linear_system_solver->b,s->linear_system_solver->n,s->linear_system_solver->x);
    
    printLinearSystem(s->linear_system_solver);     // DEBUG
    for (int i = 0; i < 4; i++)
        printf("%lf\n",s->linear_system_solver->x[i]);
    */
    
        
    return s;
}

void free_solver (struct solver_data *s)
{
    if (s->problem)
        free_problem(s->problem);

    if (s->linear_system_solver)
        free_linear_system(s->linear_system_solver);
    
    free(s);
    // TO DO
    // free_nonlinear_system
}

void Usage (int argc, char *argv[])
{
    fprintf(stderr,"===================================================================================================\n");
    fprintf(stderr,"Usage:> %s <problem_id> <linear_system_method_id> <nonlinear_system_method_id>\n",argv[0]);
    fprintf(stderr,"===================================================================================================\n");
    fprintf(stderr,"<problem_id> = Problem identifier (from Helio's handout)\n");
    fprintf(stderr,"\t1 = Problem 1\n");
    fprintf(stderr,"\t2 = Problem 2\n");
    fprintf(stderr,"\t3 = Problem 3\n");
    fprintf(stderr,"\t4 = Problem 4\n");
    fprintf(stderr,"\t5 = Problem 5\n");
    fprintf(stderr,"\t6 = Problem 6\n");
    fprintf(stderr,"\t7 = Problem 7\n");
    fprintf(stderr,"\t8 = Problem 8\n");
    fprintf(stderr,"\t9 = Problem 9\n");
    fprintf(stderr,"\t10 = Problem 10\n");
    fprintf(stderr,"<linear_system_method_id> = Linear System Solver identifier\n");
    fprintf(stderr,"\t1 = Jacobi\n");
    fprintf(stderr,"\t2 = Gauss-Seidel\n");
    fprintf(stderr,"\t3 = Conjugate Gradient\n");
    fprintf(stderr,"\t4 = BiConjugate Gradient\n");
    fprintf(stderr,"<nonlinear_system_method_id> = Nonlinear System Solver identifier\n");
    fprintf(stderr,"\t1 = Newton's method using Finite Differences\n");
    fprintf(stderr,"\t2 = Newton's method using Analitical Derivative\n");
    fprintf(stderr,"\t3 = Using Jacobi iteration\n");
    fprintf(stderr,"\t4 = Using Gauss-Seidel iteration\n");
    fprintf(stderr,"===================================================================================================\n");
}