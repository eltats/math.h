#include "s21_math.h"

long double s21_sqrt(double x) {
    if (x < 0 || x != x)
        return S21_NAN;
    if (is_inf(x) || x == 0)
        return x;
    long double start = 0, end = s21_fmax(1, x), mid = 0;
        while (s21_fabs(start - end) > 0.000001) {
            mid = (start + end) / 2;
            if (mid * mid < x)
                start = mid;
            else if (mid * mid >= x)
                end = mid;
        }
    return mid;
}
