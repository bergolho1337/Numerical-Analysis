// ************************************************************************************
// Program that solves the following 1d Poisson Equation by using the Thomas algorithm
// for tridiagonal matrix.
//
//  u_xx = f, where f = 1 and x in [0,1]
//
//  Boundary Conditions:
//
//  u(0) = 1, u(1) = 2
// ************************************************************************************


#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>

using namespace std;

#define PRINT_LINE "============================================================================================="
#define L 1.0

// Tridiagonal matrix storage
class Matrix
{
public:
    uint32_t n;     // Matrix order

    double *a;      // Sub-diagonal elements
    double *b;      // Diagonal elements
    double *c;      // Upper-diagonal elements

    double *c_star; // Auxiliary array
public:
    Matrix (const double dx);
    //~Matrix ();

    void build_finite_difference_matrix ();
    void print ();

};

// Right-hand side
class RHS 
{
public:
    uint32_t n;     // System size

    double *d;      // Elements

    double *d_star; // Auxiliary array
public:
    RHS (const double dx);

    void build_rhs (const double dx);
    void set_boundary_conditions();
};

Matrix::Matrix (const double dx)
{
    this->n = nearbyint(L/dx) + 1;

    this->a = new double[this->n];
    this->b = new double[this->n];
    this->c = new double[this->n];

    this->c_star = new double[this->n];
}

void Matrix::build_finite_difference_matrix ()
{
    for (uint32_t i = 0; i < this->n; i++)
        this->b[i] = 2.0;
    for (uint32_t i = 0; i < this->n; i++)
    {
        if (i == 0)
            this->a[i] = 0.0;
        else
            this->a[i] = -1.0;
        
        if (i == this->n-1)
            this->c[i] = 0.0;
        else
            this->c[i] = -1.0;
    }
}

void Matrix::print ()
{
    for (uint32_t i = 1; i < this->n-1; i++)
    {
        cout << a[i-1] << " " << b[i] << " " << c[i-1] << endl;
    }
}

RHS::RHS (const double dx)
{
    this->n = nearbyint(L/dx) + 1;

    this->d = new double[this->n];

    this->d_star = new double[this->n];
}

void RHS::build_rhs (const double dx)
{
    for (uint32_t i = 0; i < this->n; i++)
        this->d[i] = 1.0 * dx * dx;
}

void RHS::set_boundary_conditions ()
{
    this->d[0] = 1.0;
    this->d[this->n-1] = 2.0;
}

// Thomas algorithm for solving tridiagonal lienar system
double* solve_linear_system_using_thomas (Matrix *matrix, RHS *rhs)
{
    uint32_t n = matrix->n;
    double *a = matrix->a;
    double *b = matrix->b;
    double *c = matrix->c;
    double *c_star = matrix->c_star;
    double *d = rhs->d;
    double *d_star = rhs->d_star;

    // Output
    double *x = new double[n];

    n--; // since we start from 0 not 1
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < n; i++) 
    {
        c_star[i] = c[i] / (b[i] - a[i]*c_star[i-1]);
        d_star[i] = (d[i] - a[i]*d_star[i-1]) / (b[i] - a[i]*c_star[i-1]);
    }

    x[n] = (d[n] - a[n]*d_star[n-1]) / (b[n] - a[n]*c_star[n-1]);

    for (int i = n; i-- > 0;) 
    {
        x[i] = d_star[i] - c_star[i]*x[i+1];
    }

    return x;
}

double f (const double x)
{
    return -0.5 * powf(x,2) + 1.5*x + 1;
}

void display_usage (const char pname[])
{
    cerr << PRINT_LINE << endl;
    cerr << "Usage:> " << pname << " <dx>" << endl;
    cerr << PRINT_LINE << endl;
    cerr << "Example:> ./bin/PoissonEquation1d 0.005 > solution.dat" << endl;
    cerr << PRINT_LINE << endl;
}

int main (int argc, char **argv) 
{
    if (argc-1 != 1)
    {
        display_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    double dx = atof(argv[1]);

    struct Matrix *A = new Matrix(dx);
    struct RHS *b = new RHS(dx);

    A->build_finite_difference_matrix();
    //A->print();

    b->build_rhs(dx);
    b->set_boundary_conditions();

    double *u = solve_linear_system_using_thomas(A,b);
    for (uint32_t i = 0; i < A->n; i++)
    {
        double x = i*dx;

        cout << x << " " << u[i] << " " << f(x) << endl;
    }
  
    return 0;
}








