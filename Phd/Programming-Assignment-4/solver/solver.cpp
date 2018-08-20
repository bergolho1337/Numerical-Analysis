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
    std::cout << "\t0 = Default (1 polynomium for the whole interval)" << std::endl;
    std::cout << "\t1 = Prescribed" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "<points_file> = File with points to be interpolated" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
}