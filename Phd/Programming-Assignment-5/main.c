// ---------------------------------------------------------------------------
// Program that integrate a set of points based on a polynomial interpolation
// Created by: Lucas Berg
// ---------------------------------------------------------------------------
// How to build and run this project:
//
//  Pre-Requisites:
//      1) UNIX like system;
//      2) CMake
//      3) C Compiler (gcc or clang)
//
// 1) Execute the compilation script 
//      $ ./recompile_project.sh
// 2) Execute the program with
//      $ ./bin/Assignment5 
// ---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include "solver/solver.h"

int main (int argc, char *argv[])
{
    if (argc-1 != 4)
    {
        Usage(argc,argv);
        exit(EXIT_FAILURE);
    }

    struct solver_data *solver = new_solver_data(argc,argv);

    integrate(solver);

    free_solver_data(solver);

    return 0;
}