#include "s21_math.h"


double convertTrigonometrik(double x) {
    double res;
    if (x > 0) {
        res = x - ((int)(x / (2 * S21_M_PI)) * (2 * S21_M_PI));
    } else {
        res = x - ((int)(x / (-2 * S21_M_PI)) * (-2 * S21_M_PI));
    }
    return res;
}

long double s21_factorial(int n) {
    long double factorial = 1;
    if (n < 0) {
        factorial = 0;
    } else {
        for (; n > 1; n--) {
            factorial *= n;
        }
    }
    return factorial;
}

long double s21_sin(double x) {
    if (x == 0)
        return x;
    if (is_inf(x) || x != x)
        return S21_NAN;
    long double sin_rez = 0;
    int n = 0;
    while (x < -2*S21_M_PI || 2*S21_M_PI < x) {
        if (x > 2 * S21_M_PI) {
            x -= 2 * S21_M_PI;
        } else {
            x += 2 * S21_M_PI;
        }
    }
    while (n < 100) {
        sin_rez += s21_pow(-1, n) * s21_pow(x, 2*n+1)/s21_factorial(2*n+1);

        n++;
    }
    return sin_rez;
}

long double s21_cos(double x) {
    long double res = 0;
    if (x == 0) {
        res = 1.0;
    } else if (is_inf(x) || x != x) {
        res = S21_NAN;
    } else {
        double convert_x = convertTrigonometrik(x);
        double y = 0.00000000000000000001 + 1;
        int n = 0;
        int firstPass = 1;
        while (s21_fabs(y) > 0.00000000000000000001) {
            if (firstPass == 1) {
                y = 1;
            } else {
                y = - y * convert_x / (2 * n - 1) * convert_x / (2 * n);
            }
            firstPass = 0;
            res += y;
            n++;
        }
    }

    return res;
}

long double s21_atan(double x) {
    const long double s21_atan_1 = 0.7853981633974480L;
    if (x == 0)
        return x;
    if (x == 1)
        return s21_atan_1;
    if (x == -1)
        return -s21_atan_1;
    if (is_inf(x) && x > 0)
        return S21_M_PI_2;
    if (is_inf(x) && x < 0)
        return -S21_M_PI_2;
    if (x != x)
        return x;
    int neg = 0, per = 0, sp = 0;
    x < 0.l ? x = -x, neg = 1.l : 0;
    x > 1.l ? x = (1.l / x), per = 1 : 0;
    while (x > S21_M_PI / 50.l) {
        sp++;
        x = ((x * s21_sqrt(3)) - 1.l) * (1.l / (x + s21_sqrt(3)));
    }
    x = x * ((0.55913709l / (1.4087812l + x * x)) + 0.60310579l - 0.05160454l * (x * x));
    while (sp > 0) {
        x += (S21_P_6);
        sp--;
    }
    x = per ? (S21_M_PI_2) - x : x;
    x = neg ? -x : x;
    return x;
}

long double s21_tan(double x) {
    long double tan_rez = 0;
    if (x == 0)
        return x;
    if (is_inf(x) || x != x)
        return S21_NAN;
    tan_rez = s21_sin(x)/s21_cos(x);
    return tan_rez;
}

long double s21_asin(double x) {
    long double res = 0;
    if (x > 1.0 || x < -1.0) {
        res = S21_NAN;
    } else if (x == -1) {
        res = -S21_M_PI / 2;
    } else if (x == 1) {
        res = S21_M_PI / 2;
    } else {
        long double y = 1 + 0.00000000000000000001;
        int n = 0;
        while (s21_fabs(y) > 0.00000000000000000001) {
            if (n == 0) {
                y = x;
            } else {
                y = y * (2 * n - 1) / (2 * n + 1) * (2 * n - 1) / (2 * n) * x * x;
            }
            res += y;
            n++;
        }
    }
    return res;
}

long double s21_acos(double x) {
    long double acos_rez = 0;
    if (x == 1)
        return 0;
    acos_rez = S21_M_PI_2 - s21_asin(x);
    return acos_rez;
}
