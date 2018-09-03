#ifndef PROBLEM12_H
#define PROBLEM12_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Problem that solves the Nonlinear-System of the Assignment about the 
// Gaussian Quadruture rule passed by Heder.

// My notation for the problem ...
// w0 = x[0]
// w1 = x[1]
// w2 = x[2]
// w3 = x[3]
// t1 = x[4]
// t2 = x[5]

double f1 (const double x[]);
double f2 (const double x[]);
double f3 (const double x[]);
double f4 (const double x[]);
double f5 (const double x[]);
double f6 (const double x[]);
unsigned int getNumberEquations ();

#endif