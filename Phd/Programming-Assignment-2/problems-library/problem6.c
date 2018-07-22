#include "problem6.h"

double f1 (const double x[])
{
    return (3.0f*x[0]) - (cos(x[1]*x[2])) - (0.5f);
}

double f2 (const double x[])
{
    return (x[0]*x[0]) - (81.0f*pow(x[1]+0.1f,2)) + (sin(x[2])) + (1.06f);
}

double f3 (const double x[])
{
    return (exp(-x[0]*x[1])) + (20.0f*x[2]) + ( (10.0f*M_PI - 3.0f)/(3.0f) );
}

unsigned int getNumberEquations ()
{
    return 3;
}