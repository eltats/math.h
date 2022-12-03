#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#include <float.h>

#define S21_NAN -(0.0/0.0)
#define S21_NAN_NEGATIVE 0.0 / 0.0
#define S21_INF 1.0/0.0
#define S21_INFINITY_NEGATIVE -1.0/0.0
#define S21_INFINITY_POSITIVE 1.0/0.0
#define S21_E 2.71828182845904523536028747135266250
#define S21_EPS 1e-9
#define S21_M_PI 3.14159265358979323846264338327950288
#define S21_M_PI_2 1.57079632679489661923132169163975144
#define S21_P_6 0.52359877559
#define is_inf(x) __builtin_isinf(x)

long double s21_log(double x);
long double s21_fabs(double x);
int s21_abs(int x);
long double s21_acos(double x);
long double s21_asin(double x);
long double s21_atan(double x);
long double s21_ceil(double x);
long double s21_factorial(int n);
long double s21_cos(double x);
long double s21_exp(double x);
long double s21_floor(double x);
long double s21_fmod(double x, double y);
long double s21_pow(double base, double exp);
long double s21_sin(double x);
long double s21_sqrt(double x);
long double s21_tan(double x);
long double s21_fmax(double x, double y);

#endif  // SRC_S21_MATH_H_
