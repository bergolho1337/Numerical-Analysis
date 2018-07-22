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
    

    free_solver(solver);

    return 0;
}