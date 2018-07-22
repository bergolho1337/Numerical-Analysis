#include "problem2.h"

double f1 (const double x[])
{
    return (x[0]) + (x[0]*x[0]) - (2.0f*x[1]*x[2]) - 0.1f;
}

double f2 (const double x[])
{
    return (x[1]) - (x[1]*x[1]) + (3.0f*x[0]*x[2]) + 0.2f;
}

double f3 (const double x[])
{
    return (x[2]) + (x[2]*x[2]) + (2.0f*x[0]*x[1]) - 0.3f;
}

unsigned int getNumberEquations ()
{
    return 3;
}