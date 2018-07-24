#include "problem.h"

struct problem_data* new_problem (int problem_id)
{
    struct problem_data *p = (struct problem_data*)malloc(sizeof(struct problem_data));

    p->problem_name = getProblemName(problem_id);
    p->library_path = getLibraryPath(problem_id);
    p->neq = getNumberEquations(p->library_path);
    p->functions = getFunctions(p->handle,p->library_path,p->neq);

    fprintf(stdout,"[Problem] Solving '%s'\n",p->problem_name);
    
    return p;

}

void free_problem (struct problem_data *p)
{
    if (p->handle)
        dlclose(p->handle);

    if (p->functions)
        free(p->functions);
}

char* getProblemName (const int problem_id)
{
    switch (problem_id)
    {
        case 1: return "problem1";
        case 2: return "problem2";
        case 3: return "problem3";
        case 4: return "problem4";
        case 5: return "problem5";
        case 6: return "problem6";
        case 7: return "problem7";
        case 8: return "problem8";
        case 9: return "problem9";
        case 10: return "problem10";
    }
    fprintf(stderr,"[-] ERROR! Problem identifier not found !\n");
    return NULL;
}

char* getLibraryPath (const int problem_id)
{
    switch (problem_id)
    {
        case 1: return "./shared-libs/libproblem1.so";
        case 2: return "./shared-libs/libproblem2.so";
        case 3: return "./shared-libs/libproblem3.so";
        case 4: return "./shared-libs/libproblem4.so";
        case 5: return "./shared-libs/libproblem5.so";
        case 6: return "./shared-libs/libproblem6.so";
        case 7: return "./shared-libs/libproblem7.so";
        case 8: return "./shared-libs/libproblem8.so";
        case 9: return "./shared-libs/libproblem9.so";
        case 10: return "./shared-libs/libproblem10.so";
    }
    fprintf(stderr,"[-] ERROR! Problem identifier not found !\n");
    return NULL;
}

int getNumberEquations (const char *library_path)
{
    void *handle = dlopen(library_path,RTLD_LAZY);

    get_problem_equations_fn *neq_fn = dlsym(handle,"getNumberEquations");
    int neq = neq_fn();

    dlclose(handle);

    return neq;
}

set_problem_fn** getFunctions (void *handle, const char *library_path, const unsigned int neq)
{
    set_problem_fn **functions = (set_problem_fn**)malloc(sizeof(set_problem_fn*)*neq);

    handle = dlopen(library_path,RTLD_LAZY);
    if (!handle) 
    {
        fprintf(stderr,"%s\n",dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout,"\n[+] Problem library \"%s\" open with sucess\n",library_path);
    }
    
    char function_name[20];

    for (int i = 0; i < neq; i++)
    {
        sprintf(function_name,"f%d\0",i+1);
        functions[i] = dlsym(handle, function_name);
        if (dlerror() != NULL)  
        {
            fprintf(stderr, "\n%s function not found in the provided problem library\n",function_name);
            exit(EXIT_FAILURE);
        }
    }

    return functions;
}

void printProblem (struct problem_data *p)
{
    if (p != NULL)
    {
        double x[p->neq];
        fprintf(stdout,"Problem name = %s\n",p->problem_name);
        fprintf(stdout,"Library path = %s\n",p->library_path);
        fprintf(stdout,"Number of equations = %u\n",p->neq);
        fprintf(stdout,"Testing functions with (");
        for (int i = 0; i < p->neq-1; i++)
        {
            x[i] = 0.0f;
            fprintf(stdout,"%lf,",x[i]);
        }
        fprintf(stdout,"%lf)\n",x[p->neq-1]);
        for (int i = 0; i < p->neq; i++)
            fprintf(stdout,"Function %d = %lf\n",i+1,p->functions[i](x));
    }
    
}