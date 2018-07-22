#include "problem7.h"

double f1 (const double x[])
{
    return (x[0]) + (cos(x[0]*x[1]*x[2])) - (1.0f);
}

double f2 (const double x[])
{
    return (pow(1.0f-x[0],0.25)) + (x[1]) + (0.05f*x[2]*x[2]) - (0.15f*x[2]) - 1.0f;
}

double f3 (const double x[])
{
    return - (x[0]*x[0]) - (0.1f*x[1]*x[1]) + (0.01f*x[1]) + (x[2]) - 1.0f;
}

unsigned int getNumberEquations ()
{
    return 3;
}