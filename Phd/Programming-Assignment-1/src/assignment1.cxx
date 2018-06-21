// ***********************************************************************************
// Programming Assignment 1
//  This program calculates the value of an exponential of a value 'x' using a series
// , where we try to explore numerical problems that arise from the arithmetical 
// operations. 
//  One example is the cancelation error that occur when we try to calculate the 
// exponential of a negative number.
//
//           inf   
//           ___
//           \    (1) + (x) + (x^2 / 2!) + (x^3 / 3!) + ...
//  exp(x) = /__           
//           k=0
//
// ***********************************************************************************
 
#include <iostream>
#include <cmath>
#include <cstdio>
#include "../include/utils.h"

using namespace std;

const unsigned int NTERM = 100;      // Number of terms of the series
double factorial_table[NTERM];       // Table that stores the factorials

// Calculate the factorial table
void calcFactorialTable ()
{
    //printf("[Factorial] Calculating factorial table of %d terms ... ",NTERM);
    factorial_table[0] = 1.0;
    for (int i = 1; i < NTERM; i++)
        factorial_table[i] = i * factorial_table[i-1];
    //fflush(stdout);
    //printf("Done\n");
}

// DEBUG
void printFactorialTable ()
{
    for (int i = 0; i < NTERM; i++)
        printf("Factorial of %d = %.10lf\n",i,factorial_table[i]);
}


// Uses a naive approach to calculate the series. 
// Cancelation error will happen more often
double calcExponentialAproxNaive (const double x)
{
    //printf("[Naive] Calculating exponential aproximation of exp^(%lf) ... ",x);
    double result = 0.0;
    for (int i = 0; i < NTERM; i++)
    {
        result += pow(x,i) / factorial_table[i];
    }
    //fflush(stdout);
    //printf("Done\n");

    return result;
}

// Tries to solve the cancelation error by dividing the sum of series in two parts.
// When (x < 0), the odd terms will always be negative and the even ones positives.
// So, we sum up all the negatives and all the positives in two separate variables.
// And instead of making (NTERM / 2) operations that will suffer from cancelation errors
// we now make only 1 operation, which is the sum of positive side of the series with 
// the negative side. 
double calcExponentialAproxClever (const double x)
{
    //printf("[Clever] Calculating exponential aproximation of exp^(%lf) ... ",x);
    double sum_even = 0.0;
    double sum_odd = 0.0;
    for (int i = 0; i < NTERM; i++)
    {
        double term = pow(x,i) / factorial_table[i];
        if (i % 2 == 0)
            sum_even += term;
        else
            sum_odd += term;
    }
    //fflush(stdout);
    //printf("Done\n");
    
    return sum_even + sum_odd;
}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        printf("%s\n",PRINT_LINE);
        printf("Usage:> %s <x>\n",argv[0]);
        printf("<x> = The expoent value of the exponential to be calculated\n");
        printf("%s\n",PRINT_LINE);
        exit(EXIT_FAILURE);
    }

    // Get input from user
    double x = atof(argv[1]);

    // Calculate the factorial table
    calcFactorialTable();
    //printFactorialTable();

    // Calculate the aproximation and analitical value
    double aprox_naive = calcExponentialAproxNaive(x);
    double aprox_clever = calcExponentialAproxClever(x);
    double analit = exp(x);

    // Calculate the absolute and relative errors
    double abs_error_naive = fabs(analit-aprox_naive);
    double rel_error_naive = abs_error_naive / analit;
    double abs_error_clever = fabs(analit-aprox_clever);
    double rel_error_clever = abs_error_clever / analit;

    printInfoResults(aprox_naive,abs_error_naive,rel_error_naive,\
                    aprox_clever,abs_error_clever,rel_error_clever,\
                    analit);
    
    //writeResults(aprox_naive,abs_error_naive,rel_error_naive,\
                 aprox_clever,abs_error_clever,rel_error_clever,\
                 analit,x);
    
    return 0;
}