#ifndef SOLVER_H
#define SOLVER_H

#include "least-squares.h"
#include "linear-system.h"

static const int NEVAL = 500;

struct solver_data
{
    int n;                              // Number of phi functions to use
    int problem_id;                     // Problem identifier
    int linear_system_id;               // Linear system solver identifier
    int least_square_id;                // Least square identifier
    char *points_filename;              // Reference to the points filename
    int npoints;                        // Number of points in the input file
    double *x;                          // Coordinates from the points in the input file
    double *y;

    struct linear_system_data *linear_system;       // Pointer to the Linear system data structure
    struct least_squares_data *least_squares;       // Pointer to the Least squares data structure
};

void Usage (int argc, char *argv[]);

struct solver_data* new_solver_data (int argc, char *argv[]);
void free_solver_data (struct solver_data *s);

void solve (struct solver_data *s);
void read_points (struct solver_data *s);
double* build_matrix (struct solver_data *s);
double* build_rhs (struct solver_data *s);
double* build_matrix_2 (struct solver_data *s1);
double* build_rhs_2 (struct solver_data *s1, struct solver_data *s2);
double compute_coefficient_matrix (set_least_squares_fn *phi, const int i, const int j,\
                            const double *x, const int npoints);
double compute_coefficient_rhs (set_least_squares_fn *phi, const int i,\
                            const double *x, const double *y, const int npoints);

void solve_problem_2 (struct solver_data *s);

void print_points (struct solver_data *s);
void write_solution (struct solver_data *s);
void write_solution_problem_1 (FILE *file, struct solver_data *s);
void write_solution_problem_2 (FILE *file, struct solver_data *s);
void write_solution_problem_3 (FILE *file, struct solver_data *s);


#endif