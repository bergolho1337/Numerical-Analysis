#include "problem3.h"

double f1 (const double x[])
{
    return (10.0f*x[0]) - (2.0f*x[1]*x[1]) + (x[1]) - (2.0f*x[2]) - 5.0f;
}

double f2 (const double x[])
{
    return (8.0f*x[1]*x[1]) + (4.0f*x[2]*x[2]) - (9.0f);
}

double f3 (const double x[])
{
    return (8.0f*x[1]*x[2]) + (4.0f);
}

unsigned int getNumberEquations ()
{
    return 3;
}