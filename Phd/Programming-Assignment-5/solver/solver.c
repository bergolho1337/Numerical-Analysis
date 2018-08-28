#include "solver.h"

void Usage (int argc, char *argv[])
{
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    fprintf(stderr,"Usage:> %s <id_interpolation> <points_file> <interval_file> <id_newton_cotes>\n",argv[0]);
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    fprintf(stderr,"<id_interpolation> = Type of interpolation polynomium to use\n");
    fprintf(stderr,"\t0 = Lagrange\n");
    fprintf(stderr,"\t1 = Newton\n");
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    fprintf(stderr,"<points_file> = File with points to be integrated\n");
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    fprintf(stderr,"<interval_file> = File with interval configuration\n");
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    fprintf(stderr,"<id_newton_cotes> = Newton-Cotes integration rule to use\n");
    fprintf(stderr,"\t0 = Trapezium\n");
    fprintf(stderr,"\t1 = 1/3 Simpson\n");
    fprintf(stderr,"\t2 = 3/8 Simpson\n");
    fprintf(stderr,"\t3 = Mixed\n");
    fprintf(stderr,"** Mixed = Integration rule is choosen based on the degree of the interval **\n");
    fprintf(stderr,"** Following that we have zero error with: **\n");
    fprintf(stderr,"** Trapezium Rule -> Linear Polynomium **\n");
    fprintf(stderr,"** Simpson 1/3 Rule -> Quadratic Polynomium **\n");
    fprintf(stderr,"** Simpson 3/8 Rule -> Cubic Polynomium **\n"); 
    fprintf(stderr,"--------------------------------------------------------------------------------------------\n");
    
}

struct solver_data* new_solver_data (int argc, char *argv[])
{
    struct solver_data *s = (struct solver_data*)malloc(sizeof(struct solver_data));
    
    s->interpolation_id = atoi(argv[1]);
    s->interpolation = new_interpolation_data(s->interpolation_id,argv[2],argv[3]);

    s->newton_cotes_id = atoi(argv[4]);
    s->newton_cotes = new_newton_cotes_data(s->newton_cotes_id);

    return s;
}

void free_solver_data (struct solver_data *s)
{
    if (s->interpolation)
        free_interpolation_data(s->interpolation);

    if (s->newton_cotes)
        free_newton_cotes_data(s->newton_cotes);

    free(s);
}

void integrate (struct solver_data *s)
{
    fprintf(stdout,"\n[Solver] Calculating integral ...\n");

    // Get interpolation references
    struct interpolation_data *idata = s->interpolation;
    int ninterval = s->interpolation->ninterval;
    struct interval_data *intervals = s->interpolation->intervals;
    double *xpts = s->interpolation->x;
    double *ypts = s->interpolation->y;

    // Get Newton-Cotes references
    set_newton_cotes_fn *integrate_fn = s->newton_cotes->function;
    double total_int = 0.0;
    
    // Integrate each of the subintervals
    for (int k = 0; k < ninterval; k++)
    {
        fprintf(stdout,"--------------------------------\n");
        fprintf(stdout,"Interval %d\n",k);

        int minid = intervals[k].minid;
        int maxid = intervals[k].maxid;
        int degree = intervals[k].degree;
        double a = xpts[minid];
        double b = xpts[maxid];
        
        double local_int = integrate_fn(minid,maxid,idata,degree);

        total_int += local_int;
    }

    fprintf(stdout,"\n[Solver] Result of the integral = %.10lf\n",total_int);
}