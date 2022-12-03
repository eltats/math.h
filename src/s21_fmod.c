#include "s21_math.h"

long double s21_fmod(double x, double y) {
    if (x != x || y != y || (is_inf(x) && y == y) || (y == 0 && x == x))
        return S21_NAN;
    if (x == 0 && y != 0)
        return 0;
    if (is_inf(y) && !is_inf(x))
        return x;
    long double lx = x;
    long double ly = y;
    long long div = x / y;
    return lx - div * ly;
}
