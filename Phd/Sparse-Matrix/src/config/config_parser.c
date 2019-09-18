#include "config_parser.h"

void display_usage (char **argv) 
{

    printf ("=======================================================================================================\n");
    printf ("Usage: %s <input_matrix>\n", argv[0]);
    printf ("=======================================================================================================\n");
    printf ("Examples:\n");
    printf ("./bin/SparseMatrix inputs/Square/Non-Symmetric/sample1.mtx\n");
    printf ("=======================================================================================================\n");

    exit (EXIT_FAILURE);
}