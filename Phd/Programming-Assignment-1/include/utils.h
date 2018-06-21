#ifndef UTILS_H
#define UTILS_H

#include <cstdio>

#define PRINT_LINE "======================================================================"

void printInfoResults (const double aprox_naive, const double abs_error_naive, const double rel_error_naive,\
                       const double aprox_clever, const double abs_error_clever, const double rel_error_clever,\
                       const double analit)
{
    printf("\n");
    printf("*********** NAIVE *****************\n");
    printf("Aproximation value = %.10lf\n",aprox_naive);
    printf("Analitical value = %.10lf\n",analit);
    printf("Absolute error = %E\n",abs_error_naive);
    printf("Relative error = %E\n",rel_error_naive);
    printf("*********** NAIVE *****************\n");
    printf("\n");
    printf("*********** CLEVER *****************\n");
    printf("Aproximation value = %.10lf\n",aprox_clever);
    printf("Analitical value = %.10lf\n",analit);
    printf("Absolute error = %E\n",abs_error_clever);
    printf("Relative error = %E\n",rel_error_clever);
    printf("*********** CLEVER *****************\n");
}

void writeResults (const double aprox_naive, const double abs_error_naive, const double rel_error_naive,\
                       const double aprox_clever, const double abs_error_clever, const double rel_error_clever,\
                       const double analit, const double x)
{
    printf("%.20lf %.20lf %.20lf %.20lf %.20lf %.20lf %.20lf %.20lf\n",x,analit\
                                                      ,aprox_naive,aprox_clever\
                                                      ,abs_error_naive,abs_error_clever\
                                                      ,rel_error_naive,rel_error_clever);
}

#endif
