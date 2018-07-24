#include "nonlinear-system.h"

struct nonlinear_system_data* new_nonlinear_system (const int nonlinear_system_id)
{
    struct nonlinear_system_data *nls = (struct nonlinear_system_data*)malloc(sizeof(struct nonlinear_system_data));

    nls->method_name = getNonlinearMethodName(nonlinear_system_id);
    nls->solver = getNonlinearMethodFunction(nls->handle,nonlinear_system_id);

    return nls;
}

void free_nonlinear_system (struct nonlinear_system_data *nls)
{
    if (nls->handle)
        dlclose(nls->handle);
}

char *getNonlinearMethodName (const int nonlinear_system_id)
{
    switch (nonlinear_system_id)
    {
        case 1: return "NewtonFiniteDifference";
        case 2: return "NewtonAnaliticalDerivative";
        case 3: return "JacobiIteration";
        case 4: return "GaussSeidelIteration";
    }
}

set_nonlinear_system_fn* getNonlinearMethodFunction (void *handle, const int nonlinear_system_id)
{
    char *library_path = "./shared-libs/libdefault-nonlinear-system-solver.so";
    char *method_name = getNonlinearMethodName(nonlinear_system_id);

    handle = dlopen(library_path,RTLD_LAZY);
    if (!handle) 
    {
        fprintf(stderr,"%s\n",dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout,"\n[+] Nonlinear system library \"%s\" open with sucess\n",library_path);
    }

    set_nonlinear_system_fn *method_fn = dlsym(handle,method_name);
    if (dlerror() != NULL)  
    {
        fprintf(stderr, "[Nonlinear-System-Solver] %s function not found in the provided nonlinear system library\n",method_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout, "[Nonlinear-System-Solver] Using %s method to solve the nonlinear system\n",method_name);
    }
    return method_fn;
}
