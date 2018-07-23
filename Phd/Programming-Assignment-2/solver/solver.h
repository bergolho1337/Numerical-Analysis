#ifndef SOLVER_H
#define SOLVER_H

#include "problem.h"
#include "linear-system.h"

struct solver_data
{
    int problem_id;
    int linear_system_solver_id;
    int nonlinear_system_solver_id; 

    struct problem_data *problem;
    struct linear_system_data *linear_system_solver;
    // TODO:
    //struct nonlinear_system_data *nonlinear_system_solver;
};

struct solver_data* new_solver_data (int argc, char *argv[]);
void free_solver (struct solver_data *s);

void Usage (int argc, char *argv[]);

#endif