#ifndef SOLVER_H
#define SOLVER_H

#include "least-squares.h"
#include "linear-system.h"

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
double* build_matrix (const int n, struct least_squares_data *ls);
double* build_rhs (const int n, struct least_squares_data *ls);

void print_points (struct solver_data *s);

#endif