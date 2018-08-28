#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <dlfcn.h>

// Flag to allow a function export ...
#define EXPORT_FN

// This works like a function pointer macro ...
#define SET_INTERPOLATION(name) EXPORT_FN double name(double *x, double *y,\
                                                    const int a, const int b,\
                                                    const double z)
typedef SET_INTERPOLATION(set_interpolation_fn);

struct interval_data
{
    int degree;                     // Degree of the polynomium
    int minid;                      // Index of the first point of the interval
    int maxid;                      // Index of the last point of the interval
};

struct interpolation_data
{
    void *handle;                   // Handle to the library that interpolates the points
    char *polynomium_name;          // Name of the polynomium 
    int npoints;                    // Number of points
    int ninterval;                  // Number of intervals
    double *x;                      // X values
    double *y;                      // Y values
    struct interval_data *intervals;// Array of intervlas 
    set_interpolation_fn *function; // Pointer to the interpolation function
};

struct interpolation_data* new_interpolation_data (const int id, const char points_filename[], const char interval_filename[]);
void free_interpolation_data (struct interpolation_data *i);
void read_points (struct interpolation_data *i, const char points_filename[]);
void print_points (struct interpolation_data *i);
void read_intervals (struct interpolation_data *i, const char interval_filename[]);
void print_intervals (struct interpolation_data *i);
void read_polynomium_type (struct interpolation_data *i, const int id);
char* get_polynomium_name (const int id);
set_interpolation_fn* get_interpolation_function (void *handle, const char polynomium_name[]);


#endif



