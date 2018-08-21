// ---------------------------------------------------------------------------
// Program that interpolate a set of points
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
//      $ ./bin/Assignment4 
// ---------------------------------------------------------------------------
#include <iostream>
#include <cstdlib>
#include "solver/solver.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 != 3)
    {
        Usage(argc,argv);
        exit(EXIT_FAILURE);
    }

    Solver *solver = new Solver(argc,argv);

    solver->interpolate();

    delete solver;

    return 0;
}