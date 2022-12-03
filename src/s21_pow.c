#include "s21_math.h"

long double s21_pow(double base, double exp1) {
    long double res = 0;
    if (base < 0 && is_inf(exp1))
        return 1;
    if (base != base || exp1 != exp1)
        return S21_NAN;
    if ((base == 0 && is_inf(exp1) && exp1 < 0) \
    || (base == 0 && exp1 < 0) || (s21_fabs(base) < 1 && is_inf(exp1) && exp1 < 0) \
    || (s21_fabs(base) > 1 && is_inf(exp1) && exp1 > 0) || (is_inf(base) && exp1 > 0))
        return S21_INF;
    if (is_inf(base) && exp1 > 0 && (int)exp1 % 2 == 0)
        return -S21_INF;
    if ((base == 0 && (int)exp1 % 2 == 0 && exp1 > 0) || (s21_fabs(base) > 1 && is_inf(exp1) \
    && exp1 < 0) || (s21_fabs(base) < 1 && is_inf(exp1) && exp1 > 0) || (is_inf(base) && exp1 > 0))
        return 0;
    if (exp1 < 0) {
        base = 1/base;
        exp1 = -exp1;
    }
    if (exp1 == (int)exp1) {
        res = 1;
        for (int i = 0; i < (int)exp1; i++) {
            res *= base;
        }
        if (exp1 < 0) {
            res = 1 / res;
        }
    } else {
        res = s21_exp(exp1 * s21_log(base));
        }
    return res;
}
