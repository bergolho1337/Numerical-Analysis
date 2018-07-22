#include "problem10.h"

double f1 (const double x[])
{
    return (4.0f*x[0]) - (x[1]) + (x[2]) - (x[0]*x[3]);
}

double f2 (const double x[])
{
    return -(x[0]) + (3.0f*x[1]) - (2.0f*x[2]) - (x[1]*x[3]);
}

double f3 (const double x[])
{
    return (x[0]) - (2.0f*x[1]) + (3.0f*x[2]) - (x[2]*x[3]);
}

double f4 (const double x[])
{
    return (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]) - (1.0f);
}

unsigned int getNumberEquations ()
{
    return 4;
}