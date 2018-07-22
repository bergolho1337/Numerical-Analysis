#include "problem8.h"

double f1 (const double x[])
{
    return (x[0]) + (10.0f*x[1]);
}

double f2 (const double x[])
{
    return sqrt(5.0f)*(x[2]-x[3]);
}

double f3 (const double x[])
{
    return (pow(x[1]-x[2],2));
}

double f4 (const double x[])
{
    return sqrt(10.0f)*(pow(x[0]-x[3],2));
}

unsigned int getNumberEquations ()
{
    return 4;
}