#include "interpolation.h"

struct interpolation_data* new_interpolation_data (const int id, const char points_filename[], const char interval_filename[])
{
    struct interpolation_data *i = (struct interpolation_data*)malloc(sizeof(struct interpolation_data));

    read_points(i,points_filename);
    //print_points(i);

    read_intervals(i,interval_filename);
    //print_intervals(i);

    read_polynomium_type(i,id);
    //i->function(NULL,NULL,0,2,3.0);

    return i;
}

void free_interpolation_data (struct interpolation_data *i)
{
    if (i->handle)
        dlclose(i->handle);

    if (i->intervals)
        free(i->intervals);
    
    if (i->x)
        free(i->x);

    if (i->y)
        free(i->y);

    free(i);
}

void read_points (struct interpolation_data *i, const char points_filename[])
{
    assert(i);

    FILE *file = fopen(points_filename,"r");

    // Read number of points
    if (!fscanf(file,"%d",&i->npoints))
    {
        fprintf(stderr,"[Interpolation] Error reading points file !\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the points
    i->x = (double*)malloc(sizeof(double)*i->npoints);
    i->y = (double*)malloc(sizeof(double)*i->npoints);

    // Read values from the points
    for (int j = 0; j < i->npoints; j++)
    {
        if (!fscanf(file,"%lf %lf",&i->x[j],&i->y[j]))
        {
            fprintf(stderr,"[Interpolation] Error reading points file !\n");
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
}

void read_intervals (struct interpolation_data *i, const char interval_filename[])
{
    assert(i);

    FILE *file = fopen(interval_filename,"r");

    // Read number of intervals
    if (!fscanf(file,"%d",&i->ninterval))
    {
        fprintf(stderr,"[Interpolation] Error reading interval file !\n");
        exit(EXIT_FAILURE);
    }

    // Allocate intervals
    i->intervals = (struct interval_data*)malloc(sizeof(struct interval_data)*i->ninterval);

    // Read values from the intervals
    for (int j = 0; j < i->ninterval; j++)
    {
        if (!fscanf(file,"%d %d %d",&i->intervals[j].minid,\
                                    &i->intervals[j].maxid,\
                                    &i->intervals[j].degree))
        {
            fprintf(stderr,"[Interpolation] Error reading interval file !\n");
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);

}

void read_polynomium_type (struct interpolation_data *i, const int id)
{
    i->polynomium_name = get_polynomium_name(id);

    i->function = get_interpolation_function(i->handle,i->polynomium_name);
}

char* get_polynomium_name (const int id)
{
    switch (id)
    {
        case 0: return "Lagrange";
        case 1: return "Newton";

        default: {
                    fprintf(stderr,"[Interpolation] Error invalid interpolation method !\n");
                    exit(EXIT_FAILURE);
                }
    }
}

set_interpolation_fn* get_interpolation_function (void *handle, const char polynomium_name[])
{
    char *library_path = "./shared-libs/libdefault-interpolation-solver.so";

    handle = dlopen(library_path,RTLD_LAZY);
    if (!handle) 
    {
        fprintf(stderr,"%s\n",dlerror());
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout,"\n[Interpolation] Interpolation library \"%s\" open with sucess\n",library_path);
    }

    set_interpolation_fn *polynomium_fn = dlsym(handle,polynomium_name);
    if (dlerror() != NULL)  
    {
        fprintf(stderr, "[Interpolation] %s function not found in the provided interpolation library\n",polynomium_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stdout, "[Interpolation] Using %s method to interpolate\n",polynomium_name);
    }
    return polynomium_fn;
}

void print_points (struct interpolation_data *i)
{
    int n = i->npoints;
    double *x = i->x;
    double *y = i->y;

    fprintf(stdout,"======== POINTS ========\n");
    for (int j = 0; j < n; j++)
        fprintf(stdout,"Point %d = (%.5lf,%.5lf)\n",j,x[j],y[j]);
}

void print_intervals (struct interpolation_data *i)
{
    int n = i->ninterval;
    struct interval_data *intervals = i->intervals;

    fprintf(stdout,"======== INTERVALS ========\n");
    for (int j = 0; j < n; j++)
    {
        fprintf(stdout,"Interval %d\n",j);
        fprintf(stdout,"\tDegree = %d\n",intervals[j].degree);
        fprintf(stdout,"\tMinid = %d\n",intervals[j].minid);
        fprintf(stdout,"\tMaxid = %d\n",intervals[j].maxid);
        fprintf(stdout,"====================================\n");
    }
}