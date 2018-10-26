#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

static constexpr int neq = 2;
static constexpr double t_init = -1.0;
static constexpr double y_init[neq] = {0.0,0.0};

static constexpr double H = 1.0;
static constexpr double EPSILON = 0.05;

static constexpr int print_rate = 10;


void set_ode_initial_condition (double y[])
{
    for (int i = 0; i < neq; i++)
        y[i] = y_init[i];
}

double f1 (const double y[])
{
    return y[1];
}

double f2 (const double y[])
{
    return (H*y[0]*(1.0-y[0])*(1.0-2.0*y[0])) / (EPSILON*EPSILON);
}

double l ()
{
    return sqrt(2.0/H)*EPSILON;
}

double analit (const double t)
{
    return 0.5*(1.0 + tanh(t/(2*l()))); 
}

void write_current_solution (const double t, const double y[], FILE *out_file)
{
    fprintf(out_file,"%g ",t);
    for (int i = 0; i < neq; i++)
        fprintf(out_file,"%g ",y[i]);
    fprintf(out_file,"%g\n",analit(t));
}

void swap (double **a, double **b)
{
    double *tmp = *a;
    *a = *b;
    *b = tmp;
}

int main (int argc, char *argv[])
{
    if (argc-1 != 2)
    {
        cout << "Usage:> " << argv[0] << " <dt> <tmax>" << endl;
        cout << "<dt> = Time discretization size" << endl;
        cout << "<tmax> = Maximum time of the simulation" << endl;
        exit(EXIT_FAILURE);
    }
    double dt = atof(argv[1]);
    double tmax = atof(argv[2]);
    int M = nearbyint(tmax/dt);

    double *y_old = new double[neq];
    double *y_new = new double[neq];
    FILE *out_file = fopen("output.dat","w+");

    set_ode_initial_condition(y_old);
    for (int k = 0; k < M; k++)
    {
        double t = t_init + k*dt;

        if (k % print_rate == 0)
            write_current_solution(t,y_old,out_file);

        // Forward Euler
        y_new[0] = y_old[0] + f1(y_old)*dt;
        y_new[1] = y_old[1] + f2(y_old)*dt;

        swap(&y_old,&y_new);
    }
    fclose(out_file);
    delete [] y_old;
    delete [] y_new;

    return 0;
}