#include "lagrange.h"

double Lagrange (double *x, double *y, const int a, const int b, const int degree, const double z)
{
    //cout << "[!] Interpolating point x = " << z << " using a Lagrange polynomium of degree " << degree << endl;
    double ret = 0.0;
    for (int i = a; i <= b; i++)
    {
        double num = 1.0;
        double den = 1.0;
        for (int j = a; j <= b; j++)
        {
            if (i != j)
            {
                num *= (z - x[j]);
                den *= (x[i] - x[j]);
            }
        }
        ret += y[i] * num / den;
    }
    
    return ret;
}