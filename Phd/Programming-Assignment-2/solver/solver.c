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

    s->nonlinear_system_solver = new_nonlinear_system(s->nonlinear_system_solver_id);

    return s;
}

void free_solver (struct solver_data *s)
{
    if (s->problem)
        free_problem(s->problem);

    if (s->linear_system_solver)
        free_linear_system(s->linear_system_solver);

    if (s->nonlinear_system_solver)
        free_nonlinear_system(s->nonlinear_system_solver);
    
    free(s);
}

void solve_problem (struct solver_data *s)
{
    s->nonlinear_system_solver->x = s->nonlinear_system_solver->solver(s->linear_system_solver->A,\
                                    s->linear_system_solver->b,\
                                    s->problem->neq,\
                                    s->linear_system_solver->solver,\
                                    s->problem->functions);

    checkAnswer(s->nonlinear_system_solver->x,s->problem->neq,s->problem->functions);
}

void checkAnswer (const double *x, const int n, set_problem_fn **f)
{
    fprintf(stdout,"==============================================================================\n");
    fprintf(stdout,"[!] Checking solution for (");
    for (int i = 0; i < n-1; i++)
        fprintf(stdout,"%.10lf,",x[i]);
    fprintf(stdout,"%.10lf)\n",x[n-1]);
    fprintf(stdout,"==============================================================================\n");

    for (int i = 0; i < n; i++)
        fprintf(stdout,"Function %d = %e\n",i+1,f[i](x));

    fprintf(stdout,"==============================================================================\n");
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
    fprintf(stderr,"******* EXTRA PROBLEMS ***********\n");
    fprintf(stderr,"\t11 = Problem 11\n");
    fprintf(stderr,"\t12 = Problem 12\n");
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