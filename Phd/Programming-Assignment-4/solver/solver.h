#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>

#include "lagrange.h"

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