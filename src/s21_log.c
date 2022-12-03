#include "s21_math.h"

long double log_mercator(double x) {
    long double res = 0;
    long double y = 0.00000000000000000001 + 1;
    int n = 1;
    long double pow = 1;
    int sign = 1;
    while (s21_fabs(y) > 0.00000000000000000001) {
        pow *= x;
        y = sign * pow / n;
        res += y;
        sign *= -1;
        n++;
    }
    return res;
}

long double s21_log(double x) {
    long double res;
    if (x < 0) {
        res = S21_NAN;
    } else if (x == 0) {
        res = S21_INFINITY_NEGATIVE;
    } else if (x == S21_INFINITY_POSITIVE) {
        res = S21_INFINITY_POSITIVE;
    } else if (x != x) {
        res = x;
    } else {
        double t = (x - 1) / (x + 1);
        res = log_mercator(t) - log_mercator(-t);
    }
    return res;
}
