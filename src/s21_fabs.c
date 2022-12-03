#include "s21_math.h"

long double s21_fabs(double x) {
    long double res = 0;
    if (x != x)
        res = S21_NAN;
    else
        res = x > 0 ? x : -x;
    return res;
}
