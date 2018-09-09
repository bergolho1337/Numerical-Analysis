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
    //print_points(s);

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
    // ____________________________________________________________
    // Build the linear system structures
    s->linear_system->n = s->n;
    
    s->linear_system->A = build_matrix(s);
    //printMatrix("A",s->linear_system->A,s->linear_system->n);

    s->linear_system->b = build_rhs(s);
    //printVector("b",s->linear_system->b,s->linear_system->n);
    
    s->linear_system->x = (double*)malloc(sizeof(double)*s->n);

    // ____________________________________________________________
    // Solve the linear system
    s->linear_system->solver(s->linear_system->A,\
                             s->linear_system->b,\
                             s->linear_system->n,\
                             s->linear_system->x);

    //printVector("x",s->linear_system->x,s->linear_system->n);

    write_solution(s);
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

    // For problem 1 we need to take the natural logarithm from the y vector ... 
    if (s->problem_id == 1)
    {
        for (int i = 0; i < s->npoints; i++)
            s->y[i] = log(s->y[i]);
    }
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

double* build_matrix (struct solver_data *s)
{
    // Get the reference to the phi function ...
    int n = s->n;
    set_least_squares_fn *phi_fn = s->least_squares->function;

    // Get reference to the points dataset
    int npoints = s->npoints;
    double *x = s->x;

    // Allocate memory
    double *A = (double*)malloc(sizeof(double)*n*n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i*n+j] = compute_coefficient_matrix(phi_fn,i,j,x,npoints);
        }
    }

    return A;
}

double compute_coefficient_matrix (set_least_squares_fn *phi, const int i, const int j,\
                            const double *x, const int npoints)
{
    double ret = 0.0f;
    
    // Compute the term by: 
    //  sum_{i=0}^{npoints-1} phi_i(x_k) * phi_j(x_k)
    for (int k = 0; k < npoints; k++)
    {
        ret += phi(x[k],i)*phi(x[k],j);
    }

    return ret;
}

double* build_rhs (struct solver_data *s)
{
    // Get the reference to the phi function ...
    int n = s->n;
    set_least_squares_fn *phi_fn = s->least_squares->function;

    // Get reference to the points dataset
    int npoints = s->npoints;
    double *x = s->x;
    double *y = s->y;

    // Allocate memory
    double *b = (double*)malloc(sizeof(double)*n);

    for (int i = 0; i < n; i++)
    {
        b[i] = compute_coefficient_rhs(phi_fn,i,x,y,npoints);
    }

    return b;
}

double compute_coefficient_rhs (set_least_squares_fn *phi, const int i,\
                            const double *x, const double *y, const int npoints)
{
    double ret = 0.0f;
    
    for (int k = 0; k < npoints; k++)
    {
        ret += phi(x[k],i)*y[k];
    }   

    return ret;
}

void write_solution (struct solver_data *s)
{
    FILE *file = fopen("solution.dat","w+");

    // Output for problem 1 ...
    if (s->problem_id == 1)
    {
        double *alpha = s->linear_system->x;
        
        // First we need to recalculate the coefficient for the original problem
        double a = exp(alpha[0]);
        double b = alpha[1];
        fprintf(stdout,"\n[Solution Problem 1a] y = %.10lf * e^(%.10lf*x)\n",a,b);
        fprintf(stdout,"[Solution Problem 1b] y >= 2000 ----> x >= %.10lf\n",(log(2000.0f) - alpha[0])/alpha[1]);

        // Then, we calculate the linspace for the points to be plotted using the adjusted curve
        int npoints = s->npoints;
        double *xpts = s->x;
        //double h = (xpts[npoints-1]-xpts[0])/NEVAL;
        double h = (12.0f-xpts[0])/NEVAL;

        fprintf(file,"%d\n",NEVAL);
        // Write the adjusted curve points ...
        for (int k = 0; k < NEVAL+1; k++)
        {
            double x = xpts[0] + k*h;
            double value = a*exp(b*x);

            fprintf(file,"%.10lf %.10lf\n",x,value);
        }
    }

    fclose(file);
}