#include "newton_cotes.h"

struct newton_cotes_data* new_newton_cotes_data (const int id)
{
    struct newton_cotes_data *n = (struct newton_cotes_data*)malloc(sizeof(struct newton_cotes_data));
    
    n->integral_rule_name = get_rule_name(id);

    n->function = get_rule_function(n->handle,n->integral_rule_name);
    
    return n;
}

void free_newton_cotes_data (struct newton_cotes_data *n)
{
    if (n->handle)
        dlclose(n->handle);

    free(n);
}

char* get_rule_name (const int id)
{
    switch (id)
    {
        case 0: return "Trapezium";
        case 1: return "Simpson13";
        case 2: return "Simpson38";
        case 3: return "Mixed";
        default: {
                    fprintf(stderr,"[Newton-Cotes] Error invalid integration method !\n");
                    exit(EXIT_FAILURE);
                 }
    }
}

set_newton_cotes_fn* get_rule_function (void *handle, const char rule_name[])
{
    char *library_path = "./shared-libs/libdefault-newton-cotes-solver.so";

    handle = dlopen(library_path,RTLD_LAZY);
    if (!handle) 
    {
        fprintf(stderr,"%s\n",dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout,"\n[Newton-Cotes] Newton-Cotes library \"%s\" open with sucess\n",library_path);
    }

    set_newton_cotes_fn *newton_cotes_fn = dlsym(handle,rule_name);
    if (dlerror() != NULL)  
    {
        fprintf(stderr, "[Newton-Cotes] %s function not found in the provided Newton-Cotes library\n",rule_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout, "[Newton-Cotes] Using the %s rule to integrate the polynomium\n",rule_name);
    }
    return newton_cotes_fn;
}