#ifndef SOLVER_H
#define SOLVER_H

#include "interpolation.h"
#include "newton_cotes.h"

struct solver_data
{
    int interpolation_id;               // Interpolation identifier
    int newton_cotes_id;                // Newton-Cotes rule identifier

    struct interpolation_data *interpolation;            // Pointer to the interpolation data structure
    struct newton_cotes_data *newton_cotes;              // Pointer to the Newton-Cotes structure
};

void Usage (int argc, char *argv[]);

struct solver_data* new_solver_data (int argc, char *argv[]);
void free_solver_data (struct solver_data *s);

void integrate (struct solver_data *s);


#endif