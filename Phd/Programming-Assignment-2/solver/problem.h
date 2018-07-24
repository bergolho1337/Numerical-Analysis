#ifndef PROBLEM_H
#define PROBLEM_H

#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include <string.h>

// Flag to allow a function export ...
#define EXPORT_FN

// This works like a function pointer macro ...
#define SET_PROBLEM_TYPE(name) EXPORT_FN double name(const double x[])
typedef SET_PROBLEM_TYPE(set_problem_fn);

#define GET_PROBLEM_EQUATIONS(name) EXPORT_FN unsigned int name()
typedef GET_PROBLEM_EQUATIONS(get_problem_equations_fn);

struct problem_data
{
    void *handle;                   // Handle to the library that stores the equations of the problem
    char *problem_name;             // Problem name
    char *library_path;             // Path to the library
    unsigned int neq;               // Number of equations of the problem
    set_problem_fn **functions;     // Array of pointers to the functions(equations) of the problem
};

// Constructor and destructor
struct problem_data* new_problem (const int problem_id);
void free_problem (struct problem_data *p);

// Auxiliary functions
char* getProblemName (const int problem_id);
char* getLibraryPath (const int problem_id);
int getNumberEquations (const char *library_path);
set_problem_fn** getFunctions (void *handle, const char *library_path, const unsigned int neq);
void printProblem (struct problem_data *p);

#endif