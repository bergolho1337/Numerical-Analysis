#include "problem9.h"

double f1 (const double x[])
{
    return (x[0]) + (x[1]) - 2.0f;
}

double f2 (const double x[])
{
    return (x[0]*x[2]) + (x[1]*x[3]);
}

double f3 (const double x[])
{
    return (x[0]*x[2]*x[2]) + (x[1]*x[3]*x[3]) - (2.0f/3.0f);
}

double f4 (const double x[])
{
    return (x[0]*x[2]*x[2]*x[2]) + (x[1]*x[3]*x[3]*x[3]);
}

unsigned int getNumberEquations ()
{
    return 4;
}