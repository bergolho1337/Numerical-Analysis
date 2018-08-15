#include "problem11.h"

double f1 (const double x[])
{
    return (x[0] - 1.0f);
}

double f2 (const double x[])
{
    return (x[0]*x[1] - 1.0f);
}

unsigned int getNumberEquations ()
{
    return 2;
}