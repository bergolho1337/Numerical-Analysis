#ifndef NEWTON_COTES_H
#define NEWTON_COTES_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <dlfcn.h>

#include "interpolation.h"

// Flag to allow a function export ...
#define EXPORT_FN

// To apply the Newton-Cotes over the interpolated polynomium we need to pass 
// the indexes of the points which defines the interval plus the degree of the polynomium 
#define SET_NEWTON_COTES(name) EXPORT_FN double name(const int minid, const int maxid,\
                                                    struct interpolation_data *idata,\
                                                    const int degree)
typedef SET_NEWTON_COTES(set_newton_cotes_fn);

struct newton_cotes_data
{
    void *handle;                   // Handle to the library that integrate the points
    char *integral_rule_name;       // Name of the integral rule 
    set_newton_cotes_fn *function;  // Pointer to the Newton-Cotes function
};

struct newton_cotes_data* new_newton_cotes_data (const int id);
void free_newton_cotes_data (struct newton_cotes_data *n);
char* get_rule_name (const int id);
set_newton_cotes_fn* get_rule_function (void *handle, const char rule_name[]);


#endif