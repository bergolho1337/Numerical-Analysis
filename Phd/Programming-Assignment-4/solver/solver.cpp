#include "solver.h"

void Usage (int argc, char *argv[])
{
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "Usage:> " << argv[0] << " <id_poly> <id_interval> <points_file>" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "<id_poly> = Type of the polynomium" << std::endl;
    std::cout << "\t0 = Lagrange" << std::endl;
    std::cout << "\t1 = Newton" << std::endl;
    std::cout << "\t2 = Cubic Splines" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "<id_interval> = Type of the interval" << std::endl;
    std::cout << "\t0 = Default (the set of points is consider a single interval)" << std::endl;
    std::cout << "\t1 = Prescribed" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "<points_file> = File with points to be interpolated" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "\tLagrange polynomials with default interval:" << std::endl;
    std::cout << "$ " << argv[0] << " 0 0 ../inputs/points_1.txt\n" << std::endl;
    std::cout << "\tLagrange polynomials with prescribed interval:" << std::endl;
    std::cout << "$ " << argv[0] << " 0 1 ../inputs/points_2.txt" << std::endl;
    std::cout << "\t\t2" << std::endl;
    std::cout << "\t\t0 2" << std::endl;
    std::cout << "\t\t2 4\n" << std::endl;
    std::cout << "\tCubic splines with default interval:" << std::endl;
    std::cout << "$ " << argv[0] << " 2 0 ../inputs/points_2.txt\n" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    
}

Solver::Solver (int argc, char *argv[])
{
    polynomium_type = atoi(argv[1]);
    interval_type = atoi(argv[2]);
    points_filepath = argv[3];

    read_points();
    //print_points();

    // TO DO: Implement the CSplines using prescribed points and eliminate this if statement
    if (polynomium_type == 2 && interval_type == 1)
    {
        print_splines_message("Prescribed interval is not implemented for this method setting back to default configuration");
        interval_type = 0;
    }
    read_intervals();
    //print_intervals();

    interpolation_solver = get_interpolation_function(polynomium_type);

}

Solver::~Solver ()
{
    cout << "[!] Cleaning Solver ..." << endl;
    if (polynomium_type == 2)
        cleanup_second_derivative();
}

void Solver::interpolate ()
{
    switch (polynomium_type)
    {
        case 0: {
                    cout << "[!] Using Lagrange polynomials ..." << endl;
                    break;
                }
        case 1: {
                    cout << "[!] Using Newton polynomials ..." << endl;
                    break;
                }
        case 2: {
                    cout << "[!] Using Cubic Splines polynomials ..." << endl;
                    calc_second_derivative_csplines();
                    break;
                }
    }
    interpolate_default();
}

void Solver::interpolate_default ()
{
    // TO DO: Avoid doing this copy ...
    int n = points.size();
    double x[n], y[n];
    for (int i = 0; i < n; i++)
    {
        x[i] = points[i].getX();
        y[i] = points[i].getY();
    }

    // Interpolate the set of points
    ofstream out("scripts/interpolate_points.txt");
    int ninterval = intervals.size();
    for (int i = 0; i < ninterval; i++)
    {
        cout << "***********************************" << endl;
        cout << "Interval " << i+1 << endl;
        int minid = intervals[i].getMinId();
        int maxid = intervals[i].getMaxId(); 
        int degree = intervals[i].getDegree();
        double a = x[minid];
        double b = x[maxid];
        double h = (b - a) / NEVAL;

        // Pass through the interval
        for (int k = 0; k < NEVAL+1; k++)
        {
            double z = a + k*h;
            double value = interpolation_solver(x,y,minid,maxid,degree,z);
            out << z << " " << value << endl;
        }
    }
    
    out.close();
}

void Solver::calc_second_derivative_csplines ()
{

    // TO DO: Avoid doing this copy ...
    int n = points.size();
    double x[n], y[n];
    for (int i = 0; i < n; i++)
    {
        x[i] = points[i].getX();
        y[i] = points[i].getY();
    }

    calc_natural_splines(x,y,n);
}

void Solver::read_points ()
{
    cout << "[!] Reading points ..." << endl;

    ifstream in(points_filepath.c_str());
    uint32_t n;
    double x, y;
    in >> n;
    while (in >> x >> y)
    {
        Point p(x,y);

        points.push_back(p);
    }
    in.close();
}

void Solver::read_intervals ()
{
    switch (interval_type)
    {
        case 0: {
                    cout << "[!] Using a default interval ..." << endl;

                    Interval i(points.size(),0,points.size()-1);
                    intervals.push_back(i);
                     
                    break;
                }
        case 1: {
                    cout << "[!] Using a prescribed interval ..." << endl;

                    uint32_t ninterval, degree;
                    uint32_t a, b;

                    cout << "How many intervals you want ?" << endl;
                    cin >> ninterval;

                    for (int k = 0; k < ninterval; k++)
                    {
                        cout << "****************************************" << endl;
                        cout << "Interval " << k+1 << endl;
                        cout << "Limits [a,b]" << endl;
                        cout << "Enter the point id for a and b: ";
                        cin >> a >> b;
                        cout << "Degree of the polynomium: ";
                        cin >> degree;

                        Interval i(degree,a,b);
                        intervals.push_back(i);
                    }
                    break;
                }
    }
}

void Solver::print_intervals ()
{
    cout << "[!] Printing intervals ..." << endl;

    for (size_t i = 0; i < intervals.size(); i++)
    {
        cout << "**************************" << endl;
        cout << "Interval " << i+1 << endl;
        intervals[i].print();
        cout << "**************************" << endl;
    }
}

void Solver::print_points ()
{
    cout << "[!] Printing points ..." << endl;

    for (size_t i = 0; i < points.size(); i++)
    {
        points[i].print();
    }
}