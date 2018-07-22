#include "problem4.h"

double f1 (const double x[])
{
    return (x[0]*x[0]) + (x[1]) - 37.0f;
}

double f2 (const double x[])
{
    return (x[0]) - (x[1]*x[1]) - (5.0f);
}

double f3 (const double x[])
{
    return (x[0]) + (x[1]) + (x[2]) - (3.0f);
}

unsigned int getNumberEquations ()
{
    return 3;
}