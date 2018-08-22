#ifndef CSPLINES_H
#define CSPLINES_H

#include <iostream>
#include <cstdlib>

using namespace std;

double CSplines (double *x, double *y, const int a, const int b, const int degree, const double z);
void calc_natural_splines (const double *x, const double *y, const int n);
void cleanup_second_derivative ();
void print_splines_message (const char msg[]);

#endif