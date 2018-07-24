#ifndef SOLVER_H
#define SOLVER_H

#include "problem.h"
#include "linear-system.h"
#include "nonlinear-system.h"

struct solver_data
{
    int problem_id;                     // Problem identifier
    int linear_system_solver_id;        // Linear System method identifier
    int nonlinear_system_solver_id;     // Nonlinear System method identifier

    struct problem_data *problem;                               // Pointer to the problem data structure
    struct linear_system_data *linear_system_solver;            // Pointer to the linear system structure
    struct nonlinear_system_data *nonlinear_system_solver;      // Pointer to the nonlinear system structure
};

// Constructor and destructor
struct solver_data* new_solver_data (int argc, char *argv[]);
void free_solver (struct solver_data *s);

void solve_problem (struct solver_data *s);

// Helper functions
void Usage (int argc, char *argv[]);

#endif