#include "problem1.h"

double f1 (const double x[])
{
    return (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]) - 1.0f;
}

double f2 (const double x[])
{
    return (2.0f*x[0]*x[0]) + (x[1]*x[1]) - (4.0f*x[2]);
}

double f3 (const double x[])
{
    return (3.0f*x[0]*x[0]) - (4.0f*x[1]) + (x[2]*x[2]);
}

unsigned int getNumberEquations ()
{
    return 3;
}