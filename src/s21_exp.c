#include "s21_math.h"

long double s21_exp(double x) {
    long double res = 0;
    if (x != x) {
        res = S21_NAN;
    } else if (x == S21_INFINITY_POSITIVE) {
        res = S21_INFINITY_POSITIVE;
    } else if (x == S21_INFINITY_NEGATIVE) {
        res = 0;
    } else if (x == 0) {
        res = 1;
    } else {
        long double sum = 1, member = 1, eps = 0.00000000000000000001;
        int negative = 0;
        if (x < 0) {
            x *= -1;
            negative = 1;
        }
        for (int i = 1;; i++) {
            member *= x / i;
            sum += member;
            if (member <= eps || sum > DBL_MAX) {
                break;
            }
        }
        if (negative) {
            if (sum > DBL_MAX) {
                res = 0;
            } else {
                res = 1 / sum;
            }
        } else {
            if (sum > DBL_MAX) {
                res = S21_INFINITY_POSITIVE;
            } else {
                res = sum;
            }
        }
    }
    return res;
}
