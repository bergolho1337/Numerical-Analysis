#include <iostream>
#include <cstdlib>
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

    delete ls;


    return 0;
}