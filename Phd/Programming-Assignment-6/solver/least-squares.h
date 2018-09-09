#ifndef LEAST_SQUARES_H
#define LEAST_SQUARES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dlfcn.h>

// Flag to allow a function export ...
#define EXPORT_FN

// This works like a function pointer macro ...
#define SET_LEAST_SQUARE_SOLVER(name) EXPORT_FN double name(const double x, const int i)
typedef SET_LEAST_SQUARE_SOLVER(set_least_squares_fn);

struct least_squares_data
{
    void *handle;                   // Handle to the library that solves a linear system
    char *phi_name;                 // Name of the phi function
    set_least_squares_fn *function; // Pointer to the phi function
};

// Constructor and destructor
struct least_squares_data* new_least_squares_data (const int least_squares_id);
void free_least_squares (struct least_squares_data *ls);

// Auxiliary functions 
void getPhiFunction (struct least_squares_data *ls, const int least_squares_id);
char *getPhiFunctionName (const int least_squares_id);

#endif