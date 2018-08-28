#ifndef INTERPOLATION_H
#define INTERPOLATION_H

// Flag to allow a function export ...
#define EXPORT_FN

// This works like a function pointer macro ...
#define SET_INTERPOLATION(name) EXPORT_FN double name(double *x, double *y,\
                                                    const int a, const int b,\
                                                    const double z)
typedef SET_INTERPOLATION(set_interpolation_fn);


#endif



