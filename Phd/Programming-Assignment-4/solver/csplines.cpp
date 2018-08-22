#include "csplines.h"

// Pointer to the second derivative of each point
int nPoints;
double *s2;

// a = id_min || b = id_max 
// Currently is working with a default interval
double CSplines (double *x, double *y, const int a, const int b, const int degree, const double z)
{
    int inf, sup;
    double ret;

    // ___________________________________________________________________________________________
    // Does a binary search to locate the interval where the point 'z' is in
    if (z >= x[a] && z <= x[b])
    {
        inf = a;
        sup = b;
        while (sup - inf > 1)
        {
            int ind = (inf + sup) / 2;
            if (x[ind] > z)
                sup = ind;
            else
                inf = ind;
        }
    }
    else
    {
        ret = 0.0f;
        print_splines_message("Value to be evaluated is not is any interval !");
        exit(EXIT_FAILURE);
    }

    // ___________________________________________________________________________________________
    // Evaluate the spline using Horner procedure
    double h = x[sup] - x[inf];
    double A = (s2[sup] - s2[inf]) / (6.0f * h);
    double B = s2[inf] * 0.5f;
    double C = ((y[sup] - y[inf]) / h) - ((s2[sup] + (2.0f * s2[inf])) * (h / 6.0f));
    double D = y[inf];
    h = z - x[inf];
    ret = ((((A * h + B) * h) + C) * h) + D;
    
    return ret;
}

void calc_natural_splines (const double *x, const double *y, const int n)
{
    cout << "[CSplines] Calculating second derivatives ..." << endl;
    if (n < 3)
    {
        print_splines_message("Error ! There must be at least 3 points to interpolate");
        exit(EXIT_FAILURE);
    }
    // ___________________________________________________________________________________________
    // Allocate the global array that will store the second derivative
    nPoints = n;
    s2 = (double*)malloc(sizeof(double)*n);

    // ___________________________________________________________________________________________
    // Construct the tridiagonal simmetric system

    // Size of the system (remember boundaries will be zero here ...)
    int m = n - 2;      
    double Ha = x[1] - x[0];    
    double Deltaa = (y[1] - y[0]) / Ha;

    // Store the tridiagonal matrix using two vector (remember it is simmetric ...)
    double e[m], d[m];
    for (int i = 0; i < m; i++)
    {
        double Hb = x[i+2] - x[i+1];
        double Deltab = (y[i+2] - y[i+1]) / Hb;
        
        e[i] = Hb;
        d[i] = 2.0f * (Ha + Hb);

        s2[i+1] = 6.0f * (Deltab - Deltaa);

        Ha = Hb;
        Deltaa = Deltab;        
    }

    // ___________________________________________________________________________________________
    // Apply a Gaussian Elimination over a tridiagonal matrix
    for (int i = 1; i < m; i++)
    {
        double t = e[i-1] / d[i-1];
        d[i] -= t * e[i-1];
        s2[i+1] -= t * s2[i]; 
    }

    // ___________________________________________________________________________________________
    // Solve the system using Back-Substitution
    s2[m] /= d[m-1];
    for (int i = m-1; i >= 1; i--)
        s2[i] = (s2[i] - (e[i-1] * s2[i+1])) / d[i-1];     
    s2[0] = 0.0f;
    s2[n-1] = 0.0f;

}

void print_splines_message (const char msg[])
{
    cout << "[CSplines] " << msg << endl;
}

void cleanup_second_derivative ()
{
    free(s2);
}