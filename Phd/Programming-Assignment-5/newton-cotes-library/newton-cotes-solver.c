#include "newton-cotes-solver.h"

double Trapezium (const int minid, const int maxid, struct interpolation_data *idata, const int degree)
{
    // Get the points data for the interpolation
    double *xpts = idata->x;
    double *ypts = idata->y;

    // Calculate limits of the interval
    double a = xpts[minid];
    double b = xpts[maxid];
    fprintf(stdout,"[Newton-Cotes] a = %.10lf || b = %.10lf\n",a,b);

    // Get the interpolation function
    set_interpolation_fn *f = idata->function;

    // Calculate size of the subintervals of integral
    double h = (b-a) / NSUBINTERVAL;

    // Compute the integral following the rule ...
    double total_int = (f(xpts,ypts,minid,maxid,a) + f(xpts,ypts,minid,maxid,b)) / 2.0;
    for (int i = 1; i < NSUBINTERVAL-1; i++)
    {
        double x = a + i*h;
        total_int += f(xpts,ypts,minid,maxid,x);
    }
    total_int *= h;

    fprintf(stdout,"\tIntegral = %.10lf\n",total_int);
    return total_int;
}

double Simpson13 (const int minid, const int maxid, struct interpolation_data *idata, const int degree)
{
    // Get the points data for the interpolation
    double *xpts = idata->x;
    double *ypts = idata->y;

    // Calculate limits of the interval
    double a = xpts[minid];
    double b = xpts[maxid];
    fprintf(stdout,"[Newton-Cotes] a = %.10lf || b = %.10lf\n",a,b);

    // Get the interpolation function
    set_interpolation_fn *f = idata->function;

    // Calculate size of the subintervals
    double h = (b-a) / NSUBINTERVAL;

    // Compute the integral following the rule ...
    double total_int = (f(xpts,ypts,minid,maxid,a) + f(xpts,ypts,minid,maxid,b));
    for (int i = 1; i < NSUBINTERVAL-1; i++)
    {
        double x = a + i*h;
        if (i % 2 == 0)
            total_int += 2.0 * f(xpts,ypts,minid,maxid,x);
        else
            total_int += 4.0 * f(xpts,ypts,minid,maxid,x);
    }
    total_int *= (h / 3.0);

    fprintf(stdout,"\tIntegral = %.10lf\n",total_int);
    return total_int;
}

double Simpson38 (const int minid, const int maxid, struct interpolation_data *idata, const int degree)
{
    // Get the points data for the interpolation
    double *xpts = idata->x;
    double *ypts = idata->y;

    // Calculate limits of the interval
    double a = xpts[minid];
    double b = xpts[maxid];
    fprintf(stdout,"[Newton-Cotes] a = %.10lf || b = %.10lf\n",a,b);

    // Get the interpolation function
    set_interpolation_fn *f = idata->function;

    // Calculate size of the subintervals
    double h = (b-a) / NSUBINTERVAL;

    // Compute the integral following the rule ...
    double total_int = (f(xpts,ypts,minid,maxid,a) + f(xpts,ypts,minid,maxid,b));
    for (int i = 1; i < NSUBINTERVAL-1; i++)
    {
        double x = a + i*h;
        if (i % 3 == 0)
            total_int += 2.0 * f(xpts,ypts,minid,maxid,x);
        else
            total_int += 3.0 * f(xpts,ypts,minid,maxid,x);
    }
    total_int *= (3.0 * h / 8.0);

    fprintf(stdout,"\tIntegral = %.10lf\n",total_int);
    return total_int;
}

double Mixed (const int minid, const int maxid, struct interpolation_data *idata, const int degree)
{
    double total_int;
    switch (degree)
    {
        case 1: total_int = Trapezium(minid,maxid,idata,degree);
                break;
        case 2: total_int = Simpson13(minid,maxid,idata,degree);
                break;
        case 3: total_int = Simpson38(minid,maxid,idata,degree);
                break;
    }
    return total_int;
}