#include "least-squares.h"

struct least_squares_data* new_least_squares_data (const int least_squares_id)
{
    struct least_squares_data *ls = (struct least_squares_data*)malloc(sizeof(struct least_squares_data));

    //ls->phi_name = getPhiFunctionName(least_squares_id);
    //ls->function = getPhiFunction(ls->handle,least_squares_id);

    getPhiFunction(ls,least_squares_id);

    return ls;
}

void free_least_squares (struct least_squares_data *ls)
{
    if (ls->handle)
        dlclose(ls->handle);
    
    free(ls);
}

char *getPhiFunctionName (const int least_squares_id)
{
    switch (least_squares_id)
    {
        case 0: return "polynomial";
        default: 
        {
                fprintf(stderr,"[Least-Square] Error! Invalid identifier!\n");
                exit(EXIT_FAILURE);
        }
    }
}

void getPhiFunction (struct least_squares_data *ls, const int least_squares_id)
{
    char *library_path = "./shared-libs/libdefault-least-squares-solver.so";
    ls->phi_name = getPhiFunctionName(least_squares_id);

    ls->handle = dlopen(library_path,RTLD_LAZY);
    if (!ls->handle) 
    {
        fprintf(stderr,"%s\n",dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout,"\n[+] Least squares library \"%s\" open with sucess\n",library_path);
    }

    ls->function = dlsym(ls->handle,ls->phi_name);
    if (dlerror() != NULL)  
    {
        fprintf(stderr, "[Least-Squares-Solver] %s function not found in the provided linear system library\n",ls->phi_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout, "[Least-Squares-Solver] Using %ss to adjust the set of points\n",ls->phi_name);
    }
}
