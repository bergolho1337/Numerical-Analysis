#include "interpolation-solver.h"

double Lagrange (double *x, double *y, const int a, const int b, const double z)
{
    double ret = 0.0;
    for (int i = a; i <= b; i++)
    {
        double num = 1.0;
        double den = 1.0;
        for (int j = a; j <= b; j++)
        {
            if (i != j)
            {
                num *= (z - x[j]);
                den *= (x[i] - x[j]);
            }
        }
        ret += y[i] * num / den;
    }
    
    return ret;
}

double Newton (double *x, double *y, const int a, const int b, const double z)
{
    fprintf(stdout,"Oi eu sou o Newton !\n");
}