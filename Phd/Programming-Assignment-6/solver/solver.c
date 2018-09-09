#include "solver.h"

void Usage (int argc, char *argv[])
{
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    fprintf(stderr,"Usage:> %s <n> <points_file> <problem_id> <linear_system_id> <least_squares_id>\n",argv[0]);
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    fprintf(stderr,"<n> = Number of phi functions to use\n");
    fprintf(stderr,"\t   ___n__\n");
    fprintf(stderr,"\t   \\    /\n");
    fprintf(stderr,"\tg = \\       alpha_i * phi_i(x)\n");
    fprintf(stderr,"\t    /\n");
    fprintf(stderr,"\t   /____\\\n");
    fprintf(stderr,"\t    i = 0\n");
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    fprintf(stderr,"<points_file> = File with points to be adjusted\n");
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    fprintf(stderr,"<problem_id> = Identifier of the problem to solve\n");
    fprintf(stderr,"\t1 = Problem 1\n");
    fprintf(stderr,"\t2 = Problem 2\n");
    fprintf(stderr,"\t3 = Problem 3\n");
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    fprintf(stderr,"<linear_system_id> = Linear system method to use\n");
    fprintf(stderr,"\t0 = Jacobi\n");
    fprintf(stderr,"\t1 = Gauss_Seidel\n");
    fprintf(stderr,"\t2 = CG\n");
    fprintf(stderr,"\t3 = BiCG\n");
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    fprintf(stderr,"<least_squares_id> = Least square type\n");
    fprintf(stderr,"\t0 = phi_i(x) = x^i\n");
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    
}

struct solver_data* new_solver_data (int argc, char *argv[])
{
    struct solver_data *s = (struct solver_data*)malloc(sizeof(struct solver_data));
    
    s->n = atoi(argv[1]);
    s->points_filename = argv[2];
    s->problem_id = atoi(argv[3]);
    s->linear_system_id = atoi(argv[4]);
    s->least_square_id = atoi(argv[5]);

    s->linear_system = new_linear_system(s->linear_system_id);
    s->least_squares = new_least_squares_data(s->least_square_id);
    
    read_points(s);
    print_points(s);

    return s;
}

void free_solver_data (struct solver_data *s)
{

    if (s->linear_system)
        free_linear_system(s->linear_system);

    if (s->least_squares)
        free_least_squares(s->least_squares);    

    if (s->x)
        free(s->x);

    if (s->y)
        free(s->y);

    free(s);
}

void solve (struct solver_data *s)
{
    // Build the linear system structures
    s->linear_system->n = s->n;
    //s->linear_system->A = build_matrix(s->n,s->least_squares);
    //s->linear_system->b = build_rhs(s->n,s->least_squares);
    //s->linear_system->x = (double*)malloc(sizeof(double)*s->n);

}

void read_points (struct solver_data *s)
{
    char *input_filename = s->points_filename;
    FILE *file = fopen(input_filename,"r");

    // Read number of points
    if (!fscanf(file,"%d",&s->npoints))
        fprintf(stderr,"[-] ERROR! Reading input points file !\n");

    // Allocate memory for the points
    s->x = (double*)malloc(sizeof(double)*s->npoints);
    s->y = (double*)malloc(sizeof(double)*s->npoints);

    // Read them from the file ...
    for (int i = 0; i < s->npoints; i++)
    {
        if (!fscanf(file,"%lf %lf",&s->x[i],&s->y[i]))
            fprintf(stderr,"[-] ERROR! Reading input points file !\n");
    }

    fclose(file);
}

void print_points (struct solver_data *s)
{
    int npoints = s->npoints;
    double *x = s->x;
    double *y = s->y;

    fprintf(stdout,"***** POINTS *****\n");
    for (int i = 0; i < npoints; i++)
    {
        fprintf(stdout,"(%.10lf %10lf)\n",x[i],y[i]);
    }
}

double* build_matrix (const int n, struct least_squares_data *ls)
{
    // Get the reference to the phi function ...
    set_least_squares_fn *phi_fn = ls->function;

    // Allocate memory
    double *A = (double*)malloc(sizeof(double)*n*n);



}