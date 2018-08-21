#include "interpolation.h"

set_interpolation_fn* get_interpolation_function (const int id)
{
    set_interpolation_fn *ret;
    switch (id)
    {
        case 0: ret = Lagrange;
                break;
        case 1: ret = Newton;
                break;
        case 2: ret = CSplines;
                break;
    }
    return ret;
}