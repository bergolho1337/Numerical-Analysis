#include "linear-system.h"

struct linear_system_data* new_linear_system (const int linear_system_id)
{
    struct linear_system_data *ls = (struct linear_system_data*)malloc(sizeof(struct linear_system_data));

    ls->method_name = getMethodName(linear_system_id);
    ls->solver = getMethodFunction(ls->handle,linear_system_id);
    
    ls->n = 0;
    ls->A = NULL;
    ls->x = NULL;
    ls->b = NULL;

    return ls;
}

void free_linear_system (struct linear_system_data *ls)
{
    if (ls->handle)
        dlclose(ls->handle);

    //if (ls->method_name)
    //    free(ls->method_name);

    if (ls->A)
        free(ls->A);

    if (ls->b)
        free(ls->b);

    if (ls->x)
        free(ls->x);
    
    free(ls);
}

char *getMethodName (const int linear_system_id)
{
    switch (linear_system_id)
    {
        case 0: return "Jacobi";
        case 1: return "Gauss_Seidel";
        case 2: return "CG";
        case 3: return "BiCG";
        default: 
        {
                fprintf(stderr,"[Linear-System] Error! Invalid identifier!\n");
                exit(EXIT_FAILURE);
        }
    }
}

set_linear_system_fn* getMethodFunction (void *handle, const int linear_system_id)
{
    char *library_path = "./shared-libs/libdefault-linear-system-solver.so";
    char *method_name = getMethodName(linear_system_id);

    handle = dlopen(library_path,RTLD_LAZY);
    if (!handle) 
    {
        fprintf(stderr,"%s\n",dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout,"\n[+] Linear system library \"%s\" open with sucess\n",library_path);
    }

    set_linear_system_fn *method_fn = dlsym(handle,method_name);
    if (dlerror() != NULL)  
    {
        fprintf(stderr, "[Linear-System-Solver] %s function not found in the provided linear system library\n",method_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout, "[Linear-System-Solver] Using %s method to solve the linear system\n",method_name);
    }
    return method_fn;
}

void printLinearSystem (struct linear_system_data *ls)
{
    fprintf(stdout,"Method name = %s\n",ls->method_name);
    printMatrix("A",ls->A,ls->n);
    printVector("b",ls->b,ls->n);   
}

void printMatrix (const char *name, const double *A, const int n)
{
    fprintf(stdout,"%s\n",name);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            fprintf(stdout,"%.2lf ",A[i*n+j]);
        fprintf(stdout,"\n");
    }
}

void printVector (const char *name, const double *v, const int n)
{
    fprintf(stdout,"%s\n",name);
    for (int i = 0; i < n; i++)
    {
        fprintf(stdout,"%.2lf\n",v[i]);
    }
}