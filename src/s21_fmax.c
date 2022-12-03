#include "s21_math.h"

long double s21_fmax(double x, double y) {
    long double res = 0;
    if (x != x && y != y)
        res = S21_NAN;
    else if (x != x)
        res = y;
    else if (y != y)
        res = x;
    else
        res = x > y ? x : y;
    return res;
}
