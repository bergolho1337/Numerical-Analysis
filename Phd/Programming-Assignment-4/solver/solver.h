#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>

#include "interpolation.h"

using namespace std;

class Interval;
class Point;
class Config;
class Solver;

class Point
{
private:
    double x; 
    double y;
public:
    Point (double x, double y) : x(x), y(y) {}
    double getX () { return x; }
    double getY () { return y; }
    void print ()
    {
        cout << fixed << setprecision(10) << "(" << x << " , " << y << ")" << endl;
    }
};

class Interval
{
private:
    int degree;
    int min_id;
    int max_id;
public:
    Interval (int degree, int minid, int maxid) : degree(degree), min_id(minid), max_id(maxid) {}
    int getDegree () { return degree; }
    int getMinId () { return min_id; }
    int getMaxId () { return max_id; }
    void print ()
    {
        cout << "Degree = " << degree << endl;
        cout << "Min_id = " << min_id << endl;
        cout << "Max_id = " << max_id << endl;
    }
};

class Solver
{
private:
    uint32_t polynomium_type;
    uint32_t interval_type;
    string points_filepath;
    vector<Point> points;
    vector<Interval> intervals;
    set_interpolation_fn *interpolation_solver;

    void read_points ();
    void read_intervals ();
public:
    Solver (int argc, char *argv[]);
    ~Solver ();

    void interpolate ();
    void print_points ();
    void print_intervals ();
};

void Usage (int argc, char *argv[]);

#endif