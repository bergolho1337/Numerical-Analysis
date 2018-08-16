// -----------------------------------------------------------
// Program that implement and test different types of
// linear system solvers.
// -----------------------------------------------------------

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "../include/linear_system.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 != 3)
    {
        Usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    LinearSystem *ls = new LinearSystem(argc,argv);

    ls->solve();

    ls->write();

    ls->cond();

    delete ls;

    return 0;
}