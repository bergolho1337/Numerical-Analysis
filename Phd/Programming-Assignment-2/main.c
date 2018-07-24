// ---------------------------------------------------------------------------
// Program that solves a nonlinear system of equations
// Created by: Lucas Berg
// ---------------------------------------------------------------------------
// This project has 3 main libraries:
//
// 1) Problem library:
//      Stores the equations from the problems of Helio's handout.
// 2) Linear System library:
//      Stores some methods for solving linear system of equations.
// 3) NonLinear System library:
//      Stores some methods for solving nonlinear system of equations.
// 
// All these libraries are compiled and outputted as shared libraries with the
// '.so' extension. This was done to decrease the compilation time and also to 
// make the functions of each part of the project more readable and easy to 
// find and edit.
// ---------------------------------------------------------------------------
#include <stdio.h>
#include "solver/solver.h"

int main (int argc, char *argv[])
{
    if (argc-1 != 3)
    {
        Usage(argc,argv);
        exit(EXIT_FAILURE);
    }

    struct solver_data *solver = new_solver_data(argc,argv);
    
    solve_problem(solver);

    free_solver(solver);

    return 0;
}