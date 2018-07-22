#include "problem5.h"

double f1 (const double x[])
{
    return (x[0]*x[0]) + (2.0f*x[1]*x[1]) - (x[1]) - (2.0f*x[2]);
}

double f2 (const double x[])
{
    return (x[0]*x[0]) - (8.0f*x[1]*x[1]) + (10.0*x[2]);
}

double f3 (const double x[])
{
    return (x[0]*x[0]) / (7.0f*x[1]*x[2]) - 1.0f;
}

unsigned int getNumberEquations ()
{
    return 3;
}