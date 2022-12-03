#include <math.h>
#include <check.h>
#include <stdlib.h>
#include <limits.h>
#include "../s21_math.h"
#include "../s21_abs.c"
#include "../s21_ceil.c"
#include "../s21_floor.c"
#include "../s21_fabs.c"
#include "../s21_fmax.c"
#include "../s21_fmod.c"
#include "../s21_log.c"
#include "../s21_pow.c"
#include "../s21_sqrt.c"
#include "../s21_trigonometric.c"


START_TEST(s21_ceil_f) {
    ck_assert_ldouble_eq(s21_ceil(-15.01), ceil(-15.01));
    ck_assert_ldouble_eq(s21_ceil(15.01), ceil(15.01));
    ck_assert_ldouble_eq(s21_ceil(INFINITY), ceil(INFINITY));
    ck_assert_ldouble_eq(s21_ceil(-INFINITY), ceil(-INFINITY));
    ck_assert_ldouble_eq(s21_ceil(-0.12), ceil(-0.12));
    ck_assert_ldouble_infinite(s21_ceil(INFINITY));
    ck_assert_ldouble_infinite(ceil(INFINITY));
    ck_assert_ldouble_nan(s21_ceil(NAN));
    ck_assert_ldouble_nan(ceil(NAN));
    ck_assert_ldouble_eq(s21_ceil(0), ceil(0));
    ck_assert_ldouble_eq(s21_ceil(0.0), ceil(0.0));
    ck_assert_ldouble_eq(s21_ceil(21.21), ceil(21.21));
    ck_assert_ldouble_eq(s21_ceil(21.91), ceil(21.91));
    ck_assert_ldouble_eq(s21_ceil(-21.21), ceil(-21.21));
    ck_assert_ldouble_eq(s21_ceil(-21.91), ceil(-21.91));
    ck_assert_ldouble_eq(s21_ceil(DBL_MAX), ceil(DBL_MAX));
    double num = NAN;
    long double orig_res = ceil(num), our_res = s21_ceil(num);
    int suc = 0;
    if (isnan(orig_res) && isnan(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
} END_TEST

START_TEST(s21_floor_f) {
    ck_assert_ldouble_eq(s21_floor(0.0), floor(0.0));
    ck_assert_ldouble_eq(s21_floor(21.21), floor(21.21));
    ck_assert_ldouble_eq(s21_floor(21.91), floor(21.91));
    ck_assert_ldouble_eq(s21_floor(-21.21), floor(-21.21));
    ck_assert_ldouble_eq(s21_floor(-21.91), floor(-21.91));
    ck_assert_ldouble_infinite(s21_floor(INFINITY));
    ck_assert_ldouble_infinite(floor(INFINITY));
    ck_assert_ldouble_nan(s21_floor(NAN));
    ck_assert_ldouble_nan(floor(NAN));
    ck_assert_ldouble_eq(s21_floor(0), floor(0));
    ck_assert_ldouble_eq(s21_floor(-15.01), floor(-15.01));
    ck_assert_ldouble_eq(s21_floor(15.01), floor(15.01));
    ck_assert_ldouble_eq(s21_floor(INFINITY), floor(INFINITY));
} END_TEST

START_TEST(s21_tan_f) {
    ck_assert_ldouble_nan(s21_tan(-INFINITY));
    ck_assert_ldouble_nan(tanl(-INFINITY));
    ck_assert_ldouble_nan(s21_tan(INFINITY));
    ck_assert_ldouble_nan(tanl(INFINITY));
    ck_assert_ldouble_nan(s21_tan(NAN));
    ck_assert_ldouble_nan(tanl(NAN));
    ck_assert_double_eq_tol(s21_tan(0.0), tanl(0.0), 1e-06);
    ck_assert_double_eq_tol(s21_tan(M_PI / 6), tanl(M_PI / 6), 1e-06);
    ck_assert_double_eq_tol(s21_tan(M_PI / 4), tanl(M_PI / 4), 1e-06);
    ck_assert_double_eq_tol(s21_tan(-2 * M_PI), tanl(-2 * M_PI), 1e-06);
    ck_assert_double_eq_tol(s21_tan(M_PI), tanl(M_PI), 1e-06);
    ck_assert_double_eq_tol(s21_tan(M_PI / 3), tanl(M_PI / 3), 1e-06);
    ck_assert_double_eq_tol(s21_tan(174.532925199433), tanl(174.532925199433), 1e-06);
    ck_assert_double_eq_tol(s21_tan(-174.532925199433), tanl(-174.532925199433), 1e-06);
    for (double i = -100; i < 100; i += 2) {
        ck_assert_uint_eq(s21_tan(i), tan(i));
    }
    for (double i = -1; i < 1; i += 0.02) {
        ck_assert_uint_eq(s21_tan(i), tan(i));
    }
} END_TEST

START_TEST(s21_cos_f) {
    double value1 = 0;
    ck_assert_uint_eq(s21_cos(value1), cos(value1));
    double value2 = 112343;
    ck_assert_uint_eq(s21_cos(value2), cos(value2));
    double value3 = -312345;
    ck_assert_uint_eq(s21_cos(value3), cos(value3));
    for (double i = -S21_M_PI; i < S21_M_PI; i+= 0.01) {
        ck_assert_uint_eq(s21_cos(i), cos(i));
    }
    ck_assert_float_eq(s21_cos(0), cos(0));
    ck_assert_float_eq(s21_cos(-1), cos(-1));
    ck_assert_ldouble_nan(s21_cos(S21_INF));
    ck_assert_ldouble_nan(cosl(S21_INF));
    ck_assert_ldouble_nan(s21_cos(S21_INF));
    ck_assert_ldouble_nan(cosl(S21_INF));
    ck_assert_ldouble_nan(s21_cos(NAN));
    ck_assert_ldouble_nan(cosl(NAN));
    ck_assert_double_eq_tol(s21_cos(0.0), cosl(0.0), 1e-06);
    ck_assert_double_eq_tol(s21_cos(M_PI / 4), cosl(M_PI / 4), 1e-06);
    ck_assert_double_eq_tol(s21_cos(M_PI / 6), cosl(M_PI / 6), 1e-06);
    ck_assert_double_eq_tol(s21_cos(M_PI / 3), cosl(M_PI / 3), 1e-06);
    ck_assert_double_eq_tol(s21_cos(M_PI / 2), cosl(M_PI / 2), 1e-06);
    ck_assert_double_eq_tol(s21_cos(3 * M_PI / 2), cosl(3 * M_PI / 2), 1e-06);
    ck_assert_double_eq_tol(s21_cos(2 * M_PI), cosl(2 * M_PI), 1e-06);
    ck_assert_double_eq_tol(s21_cos(174.532925199433), cosl(174.532925199433), 1e-06);
} END_TEST

START_TEST(s21_sin_f) {
    ck_assert_ldouble_nan(s21_sin(-INFINITY));
    ck_assert_ldouble_nan(sinl(-INFINITY));
    ck_assert_ldouble_nan(s21_sin(INFINITY));
    ck_assert_ldouble_nan(sinl(INFINITY));
    ck_assert_ldouble_nan(s21_sin(NAN));
    ck_assert_ldouble_nan(sinl(NAN));
    ck_assert_double_eq_tol(s21_sin(0.0), sin(0.0), 1e-06);
    ck_assert_double_eq_tol(s21_sin(M_PI / 6), sin(M_PI / 6), 1e-06);
    ck_assert_double_eq_tol(s21_sin(M_PI / 4), sin(M_PI / 4), 1e-06);
    ck_assert_double_eq_tol(s21_sin(M_PI / 3), sin(M_PI / 3), 1e-06);
    ck_assert_double_eq_tol(s21_sin(M_PI / 2), sin(M_PI / 2), 1e-06);
    ck_assert_double_eq_tol(s21_sin(3 * M_PI / 2), sin(3 * M_PI / 2), 1e-06);
    ck_assert_double_eq_tol(s21_sin(2 * M_PI), sin(2 * M_PI), 1e-06);
    ck_assert_double_eq_tol(s21_sin(174.532925199433), sin(174.532925199433), 1e-06);
    ck_assert_double_eq_tol(s21_sin(-174.532925199433), sin(-174.532925199433), 1e-06);
    ck_assert_double_eq_tol(s21_sin(-S21_M_PI), sin(-S21_M_PI), 1e-06);
} END_TEST

START_TEST(s21_abs_f) {
    ck_assert_int_eq(s21_abs(21), abs(21));
    ck_assert_int_eq(s21_abs(-21), abs(-21));
    ck_assert_int_eq(s21_abs(-2147483647), abs(-2147483647));
    ck_assert_int_eq(s21_abs(2147483647), abs(2147483647));
    ck_assert_int_eq(s21_abs(-0), abs(-0));
    ck_assert_int_eq(s21_abs(+0), abs(+0));
    ck_assert_int_eq(s21_abs((int)NAN), abs((int)NAN));
    ck_assert_int_eq(s21_abs((int)INFINITY), abs((int)INFINITY));
    ck_assert_int_eq(s21_abs((int)-INFINITY), abs((int)-INFINITY));
    ck_assert_int_eq(s21_abs(0), abs(0));
} END_TEST

START_TEST(s21_atan_f) {
    double num = -0.99;
    long double orig_res = s21_atan(num), our_res = atan(num);
    int suc = 0;
    if ((fabsl(orig_res - our_res) <= 1e-6)) suc = 1;
    ck_assert_int_eq(1, suc);
    ck_assert_float_eq(-9999999999, -9999999999);
    ck_assert_ldouble_eq(s21_atan(INFINITY), atan(INFINITY));
    ck_assert_ldouble_eq(s21_atan(-INFINITY), atan(-INFINITY));
    double nam = NAN;
    long double orig_res1 = s21_atan(nam), our_res1 = atan(nam);
    int suc1 = 0;
    if (isnan(orig_res1) && isnan(our_res1)) suc1 = 1;
    ck_assert_int_eq(1, suc1);
    double value2 = 0.12;
    ck_assert_uint_eq(s21_atan(value2), atan(value2));
    double value3 = S21_INF;
    ck_assert_uint_eq(s21_atan(value3), atan(value3));
    ck_assert_ldouble_eq_tol(s21_atan(INFINITY), (M_PI / 2.0), 1e-6);
    ck_assert_ldouble_eq_tol(atan(INFINITY), (M_PI / 2.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_atan((-INFINITY)), ((M_PI) / -2.0), 1e-6);
    ck_assert_ldouble_eq_tol(atan(-INFINITY), ((M_PI) / -2.0), 1e-6);
    ck_assert_ldouble_nan(s21_atan(NAN));
    ck_assert_ldouble_nan(atan(NAN));
    ck_assert_ldouble_eq_tol(s21_atan(1.0), atan(1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_atan(-1.0), atan(-1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_atan(0.0), atan(0.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_atan(-0.0), atan(-0.0), 1e-6);
} END_TEST

START_TEST(s21_acos_f) {
    double value1 = 0.43;
    ck_assert_uint_eq(s21_acos(value1), acos(value1));
    double value2 = -1;
    ck_assert_uint_eq(s21_acos(value2), acos(value2));
    double value3 = -0.999;
    ck_assert_uint_eq(s21_acos(value3), acos(value3));
    for (double i = -10.; i < 10.; i += 1.) {
        ck_assert_uint_eq(s21_acos(i), acos(i));
    }
    for (double i = -1; i < 1; i += 0.01) {
        ck_assert_uint_eq(s21_acos(i), acos(i));
    }
    ck_assert_ldouble_nan(s21_acos(NAN));
    ck_assert_ldouble_nan(acos(NAN));
    ck_assert_ldouble_eq_tol(s21_acos(1.0), acos(1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(-1.0), acos(-1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(0.0), acos(0.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(-0.0), acos(-0.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(-1), acos(-1), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(1), acos(1), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(-sqrt(3) / 2), acos(-sqrt(3) / 2), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(-sqrt(2) / 2), acos(-sqrt(2) / 2), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(-1 / 2), acos(-1 / 2), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(1 / 2), acos(1 / 2), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(2 / 2), acos(2 / 2), 1e-6);
    ck_assert_ldouble_eq_tol(s21_acos(sqrt(2) / 2), acos(sqrt(2) / 2), 1e-6);
} END_TEST
START_TEST(s21_asin_f) {
    double value1 = 0.43;
    ck_assert_uint_eq(s21_asin(value1), asin(value1));
    double value2 = 0.12;
    ck_assert_uint_eq(s21_asin(value2), asin(value2));
    for (double i = 0.; i < 2; i += 1) {
        ck_assert_uint_eq(s21_asin(i), asin(i));
    }
    for (double i = -1; i < 5; i += 0.01) {
        ck_assert_uint_eq(s21_asin(i), asin(i));
    }
    ck_assert_ldouble_nan(s21_asin(INFINITY));
    ck_assert_ldouble_nan(asin(INFINITY));
    ck_assert_ldouble_nan(s21_asin(-INFINITY));
    ck_assert_ldouble_nan(asin(-INFINITY));
    ck_assert_ldouble_nan(s21_asin(NAN));
    ck_assert_ldouble_nan(asin(NAN));
    ck_assert_ldouble_eq_tol(s21_asin(1.0), asin(1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_asin(-1.0), asin(-1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_asin(0.0), asin(0.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_asin(-0.0), asin(-0.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_asin(0.0), asin(0.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_asin(-1), asin(-1), 1e-6);
    ck_assert_ldouble_eq_tol(s21_asin(1), asin(1), 1e-6);
    ck_assert_ldouble_eq_tol(s21_asin(-sqrt(3) / 2), asin(-sqrt(3) / 2), 1e-6);
    ck_assert_ldouble_eq_tol(s21_asin(-sqrt(2) / 2), asin(-sqrt(2) / 2), 1e-6);
    ck_assert_ldouble_eq_tol(s21_asin(-1 / 2), asin(-1 / 2), 1e-6);
    ck_assert_ldouble_eq_tol(s21_asin(1 / 2), asin(1 / 2), 1e-6);
} END_TEST

START_TEST(s21_fmod_f) {
  ck_assert_ldouble_eq(s21_fmod(2.34, 2.0), fmod(2.34, 2.0));
  ck_assert_ldouble_eq(s21_fmod(-2.34, 2.0), fmod(-2.34, 2.0));
  ck_assert_ldouble_eq(s21_fmod(2.34, -2.0), fmod(2.34, -2.0));
  ck_assert_ldouble_eq(s21_fmod(-2.34, -2.0), fmod(-2.34, -2.0));
  ck_assert_ldouble_eq(s21_fmod(21.21, 3), fmod(21.21, 3));
  ck_assert_ldouble_eq(s21_fmod(3, 21.21), fmod(3, 21.21));
  ck_assert_ldouble_eq(s21_fmod(-21.21, 3), fmod(-21.21, 3));
  ck_assert_ldouble_eq(s21_fmod(3, -21.21), fmod(3, -21.21));
  ck_assert_ldouble_eq(s21_fmod(100500, 9), fmod(100500, 9));
  ck_assert_ldouble_eq(s21_fmod(-100500, -9), fmod(-100500, -9));
  ck_assert_ldouble_eq(s21_fmod(-9, -100500), fmod(-9, -100500));
  ck_assert_ldouble_eq(s21_fmod(-9, -9), fmod(-9, -9));
  ck_assert_ldouble_eq(s21_fmod(10, 5), fmod(10, 5));
  ck_assert_ldouble_nan(s21_fmod(INFINITY, INFINITY));
  ck_assert_ldouble_nan(fmod(INFINITY, INFINITY));
  ck_assert_ldouble_nan(s21_fmod(NAN, NAN));
  ck_assert_ldouble_nan(fmod(NAN, NAN));
  ck_assert_ldouble_nan(s21_fmod(2.45, 0));
  ck_assert_ldouble_nan(fmod(2.45, 0));
  ck_assert_ldouble_nan(s21_fmod(0, 0));
  ck_assert_ldouble_nan(fmod(0, 0));
  ck_assert_ldouble_nan(s21_fmod(INFINITY, 0));
  ck_assert_ldouble_nan(fmod(INFINITY, 0));
} END_TEST

START_TEST(s21_exp_f) {
    double x1 = 0;
    ck_assert_int_eq(s21_exp(x1), exp(x1));
    double x2 = 1;
    ck_assert_int_eq(s21_exp(x2), exp(x2));
    double x3 = 2;
    ck_assert_int_eq(s21_exp(x3), exp(x3));
    double x4 = 3;
    ck_assert_int_eq(s21_exp(x4), exp(x4));
    double x5 = 1.5;
    ck_assert_int_eq(s21_exp(x5), exp(x5));
    double x6 = -1.5;
    ck_assert_int_eq(s21_exp(x6), exp(x6));
    double x7 = 0.5;
    ck_assert_int_eq(s21_exp(x7), exp(x7));
    double x8 = -0.5;
    ck_assert_int_eq(s21_exp(x8), exp(x8));
    double x9 = 715;
    ck_assert_int_eq(s21_exp(x9), exp(x9));
    double x10 = -715;
    ck_assert_int_eq(s21_exp(x10), exp(x10));
    double x11 = -0;
    ck_assert_int_eq(s21_exp(x11), exp(x11));
    double x12 = 'A';
    ck_assert_int_eq(s21_exp(x12), exp(x12));
    ck_assert_float_eq(s21_exp(2), exp(2));
    ck_assert_float_eq(s21_exp(-750), exp(-750));
    ck_assert_float_eq(s21_exp(-2), exp(-2));
    ck_assert_ldouble_nan(s21_exp(NAN));
    ck_assert_ldouble_nan(expl(NAN));
    ck_assert_ldouble_infinite(s21_exp(INFINITY));
  ck_assert_ldouble_infinite(expl(INFINITY));
    double num = 1000000;
    long double orig_res = exp(num), our_res = s21_exp(num);
    int suc = 0;
    if (isinf(orig_res) && isinf(our_res)) suc = 1;
    ck_assert_int_eq(1, suc);
} END_TEST

START_TEST(s21_log_f) {
    ck_assert_ldouble_infinite(s21_log(INFINITY));
    ck_assert_ldouble_infinite(log(INFINITY));
    ck_assert_ldouble_nan(s21_log(-INFINITY));
    ck_assert_ldouble_nan(log(-INFINITY));
    ck_assert_ldouble_nan(s21_log(NAN));
    ck_assert_ldouble_nan(log(NAN));
    ck_assert_ldouble_eq_tol(s21_log(1.0), log(1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_log(1.1), log(1.1), 1e-6);
    ck_assert_uint_eq(s21_log(12356), log(12356));
    ck_assert_uint_eq(s21_log(1.1), log(1.1));
    ck_assert_ldouble_eq_tol(s21_log(5.5), log(5.5), 1e-6);
    ck_assert_ldouble_eq_tol(s21_log(3.5), log(3.5), 1e-6);
    ck_assert_ldouble_eq_tol(s21_log(100), log(100), 1e-6);
    ck_assert_ldouble_eq_tol(s21_log(15.5), log(15.5), 1e-6);
    ck_assert_ldouble_eq_tol(s21_log(23.5), log(23.5), 1e-6);
    ck_assert_ldouble_infinite(s21_log(0));

    for (double i = -1.; i < 10; i+= 0.1) {
        ck_assert_uint_eq(s21_log(i), log(i));
    }
    for (double i = 0.; i < 2; i+= 0.01) {
        ck_assert_uint_eq(s21_log(i), log(i));
    }
} END_TEST

START_TEST(s21_sqrt_f) {
    ck_assert_ldouble_infinite(s21_sqrt(INFINITY));
    ck_assert_ldouble_infinite(sqrt(INFINITY));
    ck_assert_ldouble_nan(s21_sqrt(NAN));
    ck_assert_ldouble_nan(sqrt(NAN));
    ck_assert_ldouble_eq_tol(s21_sqrt(1.0), sqrt(1.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_sqrt(1.1), sqrt(1.1), 1e-6);
    ck_assert_ldouble_eq_tol(s21_sqrt(0.0), sqrt(0.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_sqrt(-0.0), sqrt(-0.0), 1e-6);
    ck_assert_ldouble_eq_tol(s21_sqrt(1), sqrt(1), 1e-6);
    ck_assert_ldouble_eq_tol(s21_sqrt(987.123), sqrt(987.123), 1e-6);
    for (double k = 0; k < 1; k+=0.001) {
        double a = s21_sqrt(k);
        double b = sqrt(k);
        ck_assert_double_eq_tol(a, b, 1e-6);
    }
    for (double k = 0; k < 100; k += 1) {
        double a = s21_sqrt(k);
        double b = sqrt(k);
        ck_assert_double_eq_tol(a, b, 1e-6);
    }
} END_TEST

START_TEST(s21_pow_f) {
    ck_assert_uint_eq(s21_pow(2.6, 3.45), pow(2.6, 3.45));
    ck_assert_uint_eq(s21_pow(3.0, 14.0), pow(3.0, 14.0));
    ck_assert_uint_eq(s21_pow(31.456, 4.3), pow(31.456, 4.3));
    ck_assert_uint_eq(s21_pow(31.456, 0.3), pow(31.456, 0.3));
    ck_assert_uint_eq(s21_pow(4.3, 4.3), pow(4.3, 4.3));
    ck_assert_uint_eq(s21_pow(-1234, 4.3), pow(-1234, 4.3));
    ck_assert_uint_eq(s21_pow(-1234, -4.3), pow(-1234, -4.3));
    ck_assert_uint_eq(s21_pow(1234, -4.3), pow(1234, -4.3));
    ck_assert_uint_eq(s21_pow(0, -4.3), pow(0, -4.3));
    ck_assert_uint_eq(s21_pow(1234, 0), pow(1234, 0));
    ck_assert_uint_eq(s21_pow(S21_NAN, -4.3), pow(S21_NAN, -4.3));
    ck_assert_uint_eq(s21_pow(1, 0), pow(1, 0));
    ck_assert_uint_eq(s21_pow(S21_INF, S21_INF),
        pow(S21_INF, S21_INF));
    ck_assert_uint_eq(s21_pow(1234, S21_INF), pow(1234, S21_INF));
    ck_assert_uint_eq(s21_pow(-1, S21_INF), pow(-1, S21_INF));
    ck_assert_uint_eq(s21_pow(-1, S21_INF), pow(-1, S21_INF));
    ck_assert_uint_eq(s21_pow(-25.0, 4.0), pow(-25.0, 4.0));
    ck_assert_uint_eq(s21_pow(-0.5, 4.0), pow(-0.5, 4.0));
    ck_assert_uint_eq(s21_pow(S21_INF, S21_INF),
        pow(S21_INF, S21_INF));
    ck_assert_uint_eq(s21_pow(S21_INF, S21_INF),
        pow(S21_INF, S21_INF));
    ck_assert_uint_eq(s21_pow(S21_INF, -5.6), pow(S21_INF, -5.6));
    ck_assert_uint_eq(s21_pow(S21_INF, S21_INF),
        pow(S21_INF, S21_INF));
    ck_assert_uint_eq(s21_pow(S21_INF, -9.0), pow(S21_INF, -9.0));
    ck_assert_uint_eq(s21_pow(S21_INF, 10.5), pow(S21_INF, 10.5));
} END_TEST


START_TEST(s21_fabs_f) {
    ck_assert_ldouble_eq(s21_fabs(0.000001), fabs(0.000001));
    ck_assert_ldouble_eq(s21_fabs(-21.000345), fabs(-21.000345));
    ck_assert_ldouble_eq(s21_fabs(-2147483600.543), fabs(-2147483600.543));
    ck_assert_ldouble_eq(s21_fabs(2147483600.4857), fabs(2147483600.4857));
    double a1 = -5.53151413431;
    ck_assert_ldouble_eq(s21_fabs(a1), fabs(a1));
    double a2 = NAN;
    ck_assert_ldouble_nan(s21_fabs(a2));
    ck_assert_ldouble_nan(fabs(a2));
    double a3 = S21_INF;
    ck_assert_ldouble_infinite(s21_fabs(a3));
    ck_assert_ldouble_infinite(fabs(a3));
    double a4 = -9519359135915.53151413431;
    ck_assert_ldouble_eq_tol(s21_fabs(a4), fabs(a4), 1e-6);
    ck_assert_ldouble_nan(s21_fabs(NAN));
    ck_assert_ldouble_nan(fabs(NAN));
    ck_assert_ldouble_eq(s21_fabs(-15.01), fabs(-15.01));
    ck_assert_ldouble_eq(s21_fabs(15.01), fabs(15.01));
    ck_assert_ldouble_eq(s21_fabs(INFINITY), fabs(INFINITY));
    ck_assert_ldouble_eq(s21_fabs(-INFINITY), fabs(-INFINITY));
}
END_TEST

START_TEST(test_abs) {
    int num = -3;
    ck_assert_double_eq_tol(s21_abs(num), abs(num), 1e-6);

    num = 3;
    ck_assert_double_eq_tol(s21_abs(num), abs(num), 1e-6);


    num = 13;
    ck_assert_double_eq_tol(s21_abs(num), abs(num), 1e-6);

    num = -11433;
    ck_assert_double_eq_tol(s21_abs(num), abs(num), 1e-6);

    num = 0;
    ck_assert_double_eq_tol(s21_abs(num), abs(num), 1e-6);

    num = -0;
    ck_assert_double_eq_tol(s21_abs(num), abs(num), 1e-6);

    num = 2147483647;
    ck_assert_double_eq_tol(s21_abs(num), abs(num), 1e-6);

    num = -2147483646;
    ck_assert_double_eq_tol(s21_abs(num), abs(num), 1e-6);

    num = INT_MIN;
    ck_assert_double_eq_tol(s21_abs(num), abs(num), 1e-6);

    num = INT_MAX;
    ck_assert_double_eq_tol(s21_abs(num), abs(num), 1e-6);
} END_TEST

START_TEST(test_acos_1) {
    double num = -3;
    ck_assert_double_nan(s21_acos(num));
}
START_TEST(test_acos_2) {
    double num = 3;
    ck_assert_ldouble_nan(s21_acos(num));
}
START_TEST(test_acos_3) {
    ck_assert_ldouble_nan(s21_acos(NAN));
}
START_TEST(test_acos_4) {
    double num = 0;
    ck_assert_ldouble_eq_tol(s21_acos(num), acos(num), 1e-09);
}
START_TEST(test_acos_5) {
    double num = 1;
    ck_assert_ldouble_eq_tol(s21_acos(num), acos(num), 1e-09);
}
START_TEST(test_acos_6) {
    double num = 1.01;
    ck_assert_ldouble_nan(s21_acos(num));
}
START_TEST(test_acos_7) {
    double num = -0;
    ck_assert_ldouble_eq_tol(s21_acos(num), acos(num), 1e-09);
}
START_TEST(test_acos_8) {
    double num = 0.9;
    ck_assert_ldouble_eq_tol(s21_acos(num), acos(num), 1e-09);
}
START_TEST(test_acos_9) {
    double num = 0.29;
    ck_assert_ldouble_eq_tol(s21_acos(num), acos(num), 1e-09);
}
START_TEST(test_acos_10) {
    double num = 0.34;
    ck_assert_ldouble_eq_tol(s21_acos(num), acos(num), 1e-09);
}
START_TEST(test_acos_11) {
    double num = -0.29;
    ck_assert_ldouble_eq_tol(s21_acos(num), acos(num), 1e-09);
}
START_TEST(test_acos_12) {
    double num = -0.34;
    ck_assert_ldouble_eq_tol(s21_acos(num), acos(num), 1e-09);
}
START_TEST(test_acos_13) {
    double num = -0.9;
    ck_assert_ldouble_eq_tol(s21_acos(num), acos(num), 1e-09);
}
START_TEST(test_acos_14) {
    double num = 0.49;
    ck_assert_ldouble_eq_tol(s21_acos(num), acos(num), 1e-09);
}
START_TEST(test_acos_15) {
    double num = -0.22;
    ck_assert_ldouble_eq_tol(s21_acos(num), acos(num), 1e-09);
}

START_TEST(test_asin_1) {
    double num = -3;
    ck_assert_ldouble_nan(s21_asin(num));
}
START_TEST(test_asin_2) {
    double num = 3;
    ck_assert_ldouble_nan(s21_asin(num));
}
START_TEST(test_asin_3) {
    double num = 0;
    ck_assert_ldouble_eq_tol(s21_asin(num), asin(num), 1e-6);
}
START_TEST(test_asin_4) {
    double num = 0.9;
    ck_assert_ldouble_eq_tol(s21_asin(num), asin(num), 1e-6);
}
START_TEST(test_asin_5) {
    double num = -0.9;
    ck_assert_ldouble_eq_tol(s21_asin(num), asin(num), 1e-6);
}
START_TEST(test_asin_6) {
    double num = 0.49;
    ck_assert_ldouble_eq_tol(s21_asin(num), asin(num), 1e-6);
}
START_TEST(test_asin_7) {
    double num = -0.22;
    ck_assert_ldouble_eq_tol(s21_asin(num), asin(num) , 1e-6);
}
START_TEST(test_asin_8) {
    double num = -0.31;
    ck_assert_ldouble_eq_tol(s21_asin(num), asin(num) , 1e-6);
}
START_TEST(test_asin_9) {
    double num = -0.67;
    ck_assert_ldouble_eq_tol(s21_asin(num), asin(num) , 1e-6);
}
START_TEST(test_asin_10) {
    double num = 0.53;
    ck_assert_ldouble_eq_tol(s21_asin(num), asin(num) , 1e-6);
}
START_TEST(test_asin_11) {
    ck_assert_ldouble_nan(s21_asin(NAN));
}

START_TEST(test_log_1) {
    double num = 3;
    ck_assert_double_eq_tol(s21_log(num), log(num), 1e-6);
}
START_TEST(test_log_2) {
    double num = 1;
    ck_assert_double_eq_tol(s21_log(num), log(num), 1e-6);
}
START_TEST(test_log_3) {
    double num = 0;
    ck_assert_double_infinite(s21_log(num));
}
START_TEST(test_log_4) {
    double num = -0;
    ck_assert_double_infinite(s21_log(num));
}
START_TEST(test_log_5) {
    double num = 0.99;
    ck_assert_double_eq_tol(s21_log(num), log(num), 1e-6);
}
START_TEST(test_log_6) {
    double num = -0.99;
    ck_assert_double_nan(s21_log(num));
}
START_TEST(test_log_7) {
    double num = 0.49;
    ck_assert_double_eq_tol(s21_log(num), log(num), 1e-6);
}
START_TEST(test_log_8) {
    double num = 1200405.35;
    ck_assert_double_eq_tol(s21_log(num), log(num), 1e-6);
    num = 12005.5;
    ck_assert_double_eq_tol(s21_log(num), log(num), 1e-6);
    num = 125;
    ck_assert_double_eq_tol(s21_log(num), log(num), 1e-6);
    num = 854.3;
    ck_assert_double_eq_tol(s21_log(num), log(num), 1e-6);
    num = 1214.1253;
    ck_assert_double_eq_tol(s21_log(num), log(num), 1e-6);
}
START_TEST(test_log_9) {
    double num = -131505.312;
    ck_assert_double_nan(s21_log(num));
}
START_TEST(test_log_10) {
    double num = -0.22;
    ck_assert_double_nan(s21_log(num));
}
START_TEST(test_log_11) {
    ck_assert_double_nan(s21_log(NAN));
}
START_TEST(test_log_12) {
    ck_assert_ldouble_infinite(s21_log(INFINITY));
}

START_TEST(test_tan_1) {
    ck_assert_double_eq_tol(s21_tan(-1), tan(-1), 1e-6);
    ck_assert_double_eq_tol(s21_tan(1), tan(1), 1e-6);
    ck_assert_double_eq_tol(s21_tan(10), tan(10), 1e-6);
    ck_assert_double_eq_tol(s21_tan(0), tan(0), 1e-6);
    ck_assert_double_eq_tol(s21_tan(0x14BD), tan(0x14BD), 1e-6);
    ck_assert_double_eq_tol(s21_tan(145), tan(145), 1e-6);
    ck_assert_double_eq_tol(s21_tan(16), tan(16), 1e-6);
    ck_assert_double_eq_tol(s21_tan(-16), tan(-16), 1e-6);
}
START_TEST(test_tan_2) {
    double num = 3;
    ck_assert_ldouble_eq_tol(s21_tan(num), tan(num), 1e-6);
}
START_TEST(test_tan_3) {
    double num = 0;
    ck_assert_ldouble_eq_tol(s21_tan(num), tan(num), 1e-6);
}
START_TEST(test_tan_4) {
    double num = -0;
    ck_assert_ldouble_eq_tol(s21_tan(num), tan(num), 1e-6);
}
START_TEST(test_tan_5) {
    double num = 0.9;
    ck_assert_ldouble_eq_tol(s21_tan(num), tan(num), 1e-6);
}
START_TEST(test_tan_6) {
    double num = 1;
    ck_assert_ldouble_eq_tol(s21_tan(num), tan(num), 1e-6);
}
START_TEST(test_tan_7) {
    double num = -0.9;
    ck_assert_ldouble_eq_tol(s21_tan(num), tan(num), 1e-6);
}
START_TEST(test_tan_8) {
    double num = 0.49;
    ck_assert_ldouble_eq_tol(s21_tan(num), tan(num), 1e-6);
}
START_TEST(test_tan_9) {
    double num = -0.22;
    ck_assert_ldouble_eq_tol(s21_tan(num), tan(num) , 1e-6);
}
START_TEST(test_tan_10) {
    ck_assert_ldouble_nan(s21_tan(NAN));
}
START_TEST(test_tan_11) {
    ck_assert_ldouble_nan(s21_tan(INFINITY));
}
START_TEST(test_tan_12) {
    ck_assert_ldouble_nan(s21_tan(-INFINITY));
}

START_TEST(test_atan_1) {
    double num = 0.5;
    ck_assert_ldouble_eq_tol(s21_atan(num), atan(num), 1e-6);
    ck_assert_ldouble_eq_tol(s21_atan(INFINITY), atan(INFINITY), 1e-6);
    ck_assert_ldouble_eq_tol(s21_atan(-INFINITY), atan(-INFINITY), 1e-6);
}
START_TEST(test_atan_2) {
    double num = 3;
    ck_assert_ldouble_eq_tol(s21_atan(num), atan(num), 1e-6);
}
START_TEST(test_atan_3) {
    double num = 0;
    ck_assert_ldouble_eq_tol(s21_atan(num), atan(num), 1e-6);
}
START_TEST(test_atan_4) {
    double num = -0;
    ck_assert_ldouble_eq_tol(s21_atan(num), atan(num), 1e-6);
}
START_TEST(test_atan_5) {
    double num = 1;
    ck_assert_ldouble_eq_tol(s21_atan(num), atan(num), 1e-6);
}
START_TEST(test_atan_6) {
    double num = -2;
    ck_assert_ldouble_eq_tol(s21_atan(num), atan(num), 1e-6);
}
START_TEST(test_atan_7) {
    double num = 2;
    ck_assert_ldouble_eq_tol(s21_atan(num), atan(num), 1e-6);
}
START_TEST(test_atan_8) {
    double num = 0.9;
    ck_assert_ldouble_eq_tol(s21_atan(num), atan(num), 1e-6);
}
START_TEST(test_atan_9) {
    ck_assert_ldouble_nan(s21_atan(NAN));
}
START_TEST(test_atan_10) {
    double num = 0.5;
    ck_assert_ldouble_eq_tol(s21_atan(num), atan(num) , 1e-6);
}
START_TEST(test_atan_11) {
    ck_assert_ldouble_nan(s21_atan(NAN));
}

START_TEST(test_ceil_1) {
    double num = -3;
    ck_assert_ldouble_eq_tol(s21_ceil(num), ceil(num), 1e-6);
}
START_TEST(test_ceil_2) {
    double num = 3;
    ck_assert_ldouble_eq_tol(s21_ceil(num), ceil(num), 1e-6);
}
START_TEST(test_ceil_3) {
    double num = 0;
    ck_assert_ldouble_eq_tol(s21_ceil(num), ceil(num), 1e-6);
}
START_TEST(test_ceil_4) {
    double num = -0;
    ck_assert_ldouble_eq_tol(s21_ceil(num), ceil(num), 1e-6);
}
START_TEST(test_ceil_5) {
    double num = 1;
    ck_assert_ldouble_eq_tol(s21_ceil(num), ceil(num), 1e-6);
}
START_TEST(test_ceil_6) {
    double num = -5;
    ck_assert_ldouble_eq_tol(s21_ceil(num), ceil(num), 1e-6);
}
START_TEST(test_ceil_7) {
    double num = -2;
    ck_assert_ldouble_eq_tol(s21_ceil(num), ceil(num), 1e-6);
}
START_TEST(test_ceil_8) {
    double num = -0.22;
    ck_assert_ldouble_eq_tol(s21_ceil(num), ceil(num) , 1e-6);
}
START_TEST(test_ceil_9) {
    double num = -15;
    ck_assert_ldouble_eq_tol(s21_ceil(num), ceil(num), 1e-6);
}
START_TEST(test_ceil_10) {
    double num = 15;
    ck_assert_ldouble_eq_tol(s21_ceil(num), ceil(num), 1e-6);
}
START_TEST(test_ceil_11) {
    ck_assert_ldouble_nan(s21_ceil(NAN));
}
START_TEST(test_ceil_12) {
    ck_assert_ldouble_infinite(s21_ceil(INFINITY));
}
START_TEST(test_ceil_13) {
    ck_assert_ldouble_infinite(s21_ceil(-INFINITY));
}
START_TEST(test_ceil_14) {
    double num = 6;
    ck_assert_ldouble_eq_tol(s21_ceil(num), ceil(num), 1e-6);
}
START_TEST(test_ceil_15) {
    double num = -0.99;
    ck_assert_ldouble_eq_tol(s21_ceil(num), ceil(num), 1e-6);
}

START_TEST(test_floor_1) {
    double num = -3;
    ck_assert_ldouble_eq_tol(s21_floor(num), floor(num), 1e-6);
}
START_TEST(test_floor_2) {
    double num = 3;
    ck_assert_ldouble_eq_tol(s21_floor(num), floor(num), 1e-6);
}
START_TEST(test_floor_3) {
    double num = 0;
    ck_assert_ldouble_eq_tol(s21_floor(num), floor(num), 1e-6);
}
START_TEST(test_floor_4) {
    double num = -0;
    ck_assert_ldouble_eq_tol(s21_floor(num), floor(num), 1e-6);
}
START_TEST(test_floor_5) {
    double num = 1;
    ck_assert_ldouble_eq_tol(s21_floor(num), floor(num), 1e-6);
}
START_TEST(test_floor_6) {
    double num = -5;
    ck_assert_ldouble_eq_tol(s21_floor(num), floor(num), 1e-6);
}
START_TEST(test_floor_7) {
    double num = -2;
    ck_assert_ldouble_eq_tol(s21_floor(num), floor(num), 1e-6);
}
START_TEST(test_floor_8) {
    double num = 2;
    ck_assert_ldouble_eq_tol(s21_floor(num), floor(num) , 1e-6);
}
START_TEST(test_floor_9) {
    double num = DBL_MIN;
    ck_assert_ldouble_eq_tol(s21_floor(num), floor(num), 1e-6);
}
START_TEST(test_floor_10) {
    double num = 15;
    ck_assert_ldouble_eq_tol(s21_floor(num), floor(num), 1e-6);
}
START_TEST(test_floor_11) {
    ck_assert_ldouble_nan(s21_floor(NAN));
}
START_TEST(test_floor_12) {
    ck_assert_ldouble_infinite(s21_floor(INFINITY));
}
START_TEST(test_floor_13) {
    ck_assert_ldouble_infinite(s21_floor(-INFINITY));
}
START_TEST(test_floor_14) {
    double num = 6;
    ck_assert_ldouble_eq_tol(s21_floor(num), floor(num), 1e-6);
}
START_TEST(test_floor_15) {
    double num = -0.9;
    ck_assert_ldouble_eq_tol(s21_floor(num), floor(num), 1e-6);
}

START_TEST(test_cos_3) {
    double num = 0;
    ck_assert_ldouble_eq_tol(s21_cos(num), cos(num), 1e-6);
}
START_TEST(test_cos_4) {
    double num = -0;
    ck_assert_ldouble_eq_tol(s21_cos(num), cos(num), 1e-6);
}
START_TEST(test_cos_9) {
    ck_assert_ldouble_nan(s21_cos(NAN));
}
START_TEST(test_cos_10) {
    ck_assert_ldouble_nan(s21_acos(INFINITY));
}
START_TEST(test_cos_11) {
    ck_assert_ldouble_nan(s21_acos(-INFINITY));
}

START_TEST(test_exp_1) {
    double num = -3;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_2) {
    double num = 3;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_3) {
    double num = 0;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_4) {
    double num = -0;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_5) {
    double num = 0.99;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_6) {
    double num = 1;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_7) {
    double num = -0.99;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_8) {
    double num = 1;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_9) {
    double num = 0.49;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_10) {
    double num = 6;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_11) {
    double num = -2;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_12) {
    double num = 14;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_13) {
    double num = 0.1;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_14) {
    double num = 10;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_15) {
    double num = DBL_MIN;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num), 1e-6);
}
START_TEST(test_exp_16) {
    double num = -0.22;
    ck_assert_ldouble_eq_tol(s21_exp(num), exp(num) , 1e-6);
}
START_TEST(test_exp_17) {
    ck_assert_ldouble_nan(s21_exp(NAN));
}
START_TEST(test_exp_18) {
    ck_assert_ldouble_infinite(s21_exp(INFINITY));
}
START_TEST(test_exp_19) {
    ck_assert_ldouble_eq_tol(s21_exp(-INFINITY), exp(-INFINITY), 1e-6);
}

START_TEST(test_fmod_1) {
    double num = -16;
    double y = 3;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_2) {
    double num = 3;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_3) {
    double num = 0;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_4) {
    double num = -0;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_5) {
    double num = 2147483648;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_6) {
    double num = INT_MIN;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_7) {
    double num = INT_MAX;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_8) {
    double num = INT_MAX;
    double y = INT_MAX;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_9) {
    double num = INT_MAX;
    double y = INT_MIN;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_10) {
    double num = INT_MIN;
    double y = INT_MIN;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_11) {
    double num = INT_MAX;
    double y = INT_MIN;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_12) {
    double num = 0;
    double y = 124.236;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_14) {
    double num = 77;
    double y = 0;
    ck_assert_double_nan(s21_fmod(num, y));
}
START_TEST(test_fmod_15) {
    double num = 77;
    double y = INFINITY;
    ck_assert_double_eq_tol(s21_fmod(num, y), fmod(num, y), 1e-6);
}
START_TEST(test_fmod_16) {
    double num = 77;
    double y = 0;
    ck_assert_double_nan(s21_fmod(num, y));
}
START_TEST(test_fmod_17) {
    double num = NAN;
    double y = 0;
    ck_assert_double_nan(s21_fmod(num, y));
    num = 0;
    y = NAN;
    ck_assert_double_nan(s21_fmod(num, y));
    num = 1;
    y = NAN;
    ck_assert_double_nan(s21_fmod(num, y));
    num = NAN;
    y = 1;
    ck_assert_double_nan(s21_fmod(num, y));
    num = INFINITY;
    y = 1;
    ck_assert_double_nan(s21_fmod(num, y));
}

START_TEST(test_pow_1) {
    double num = 3;
    double y = 3;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
}
START_TEST(test_pow_2) {
    double num = 3;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
    num = 2;
    y = -3;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
    num = -2;
    y = -3;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
}
START_TEST(test_pow_3) {
    double num = 100;
    double y = 5;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
}
START_TEST(test_pow_4) {
    double num = 124;
    double y = 97;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
}
START_TEST(test_pow_5) {
    double num = 3.1435;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
}
START_TEST(test_pow_6) {
    double num = 3.1435;
    double y = 8.2;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
}
START_TEST(test_pow_7) {
    double num = 0;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
}
START_TEST(test_pow_8) {
    double num = -0;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
}
START_TEST(test_pow_9) {
    double num = 2148;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
}
START_TEST(test_pow_10) {
    double num = -2147483646;
    double y = 1.2;
    ck_assert_double_nan(s21_pow(num, y));
}
START_TEST(test_pow_12) {
    double num = 10;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
}
START_TEST(test_pow_17) {
    double num = 0;
    double y = 124.236;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
}
START_TEST(test_pow_18) {
    double num = INFINITY;
    double y = 124.236;
    ck_assert_double_infinite(s21_pow(num, y));
}
START_TEST(test_pow_19) {
    double num = 77;
    double y = 0;
    ck_assert_double_eq_tol(s21_pow(num, y), pow(num, y), 1e-6);
}
START_TEST(test_pow_20) {
    double num = 77;
    double y = INFINITY;
    ck_assert_double_infinite(s21_pow(num, y));
}
START_TEST(test_pow_21) {
    double num = NAN;
    double y = 0;
    ck_assert_double_nan(s21_pow(num, y));
}



START_TEST(test_sin) {
    double num = 0;

    ck_assert_double_eq_tol(s21_sin(num), sin(num), 1e-6);

    num = -0;
    ck_assert_double_eq_tol(s21_sin(num), sin(num), 1e-6);

    ck_assert_double_nan(s21_sin(NAN));
} END_TEST

START_TEST(test_sqrt) {
    double num = -3;
    ck_assert_double_nan(s21_sqrt(num));

    num = 3;
    ck_assert_ldouble_eq_tol(s21_sqrt(num), sqrt(num), 1e-6);

    num = 0;
    ck_assert_ldouble_eq_tol(s21_sqrt(num), sqrt(num), 1e-6);

    num = -0;
    ck_assert_ldouble_eq_tol(s21_sqrt(num), sqrt(num), 1e-6);

    num = 0.99;
    ck_assert_ldouble_eq_tol(s21_sqrt(num), sqrt(num), 1e-6);

    num = -0.99;
    ck_assert_double_nan(s21_sqrt(num));

    num = 0.49;
    ck_assert_ldouble_eq_tol(s21_sqrt(num), sqrt(num), 1e-6);


    num = 214352;
    ck_assert_ldouble_eq_tol(s21_sqrt(num), sqrt(num), 1e-6);

    num = -214352;
    ck_assert_double_nan(s21_sqrt(num));

    num = 241.2146574;
    ck_assert_ldouble_eq_tol(s21_sqrt(num), sqrt(num), 1e-6);

    num = -241.2146574;
    ck_assert_double_nan(s21_sqrt(num));

    num = -0.22;
    ck_assert_double_nan(s21_sqrt(num));

    ck_assert_double_nan(s21_sqrt(NAN));
} END_TEST

START_TEST(test_fabs_1) {
    double num = 3;
    ck_assert_double_eq_tol(s21_fabs(num), fabs(num), 1e-6);
}
START_TEST(test_fabs_2) {
    double num = 0;
    ck_assert_double_eq_tol(s21_fabs(num), fabs(num), 1e-6);
}
START_TEST(test_fabs_3) {
    double num = -0;
    ck_assert_double_eq_tol(s21_fabs(num), fabs(num), 1e-6);
}
START_TEST(test_fabs_4) {
    double num = 2147483648;
    ck_assert_double_eq_tol(s21_fabs(num), fabs(num), 1e-6);
}
START_TEST(test_fabs_5) {
    double num = -2147483646;
    ck_assert_double_eq_tol(s21_fabs(num), fabs(num), 1e-6);
}
START_TEST(test_fabs_6) {
    double num = INT_MIN;
    ck_assert_double_eq_tol(s21_fabs(num), fabs(num), 1e-6);
}
START_TEST(test_fabs_7) {
    double num = INT_MAX;
    ck_assert_double_eq_tol(s21_fabs(num), fabs(num), 1e-6);
}
START_TEST(test_fabs_8) {
    ck_assert_double_infinite(s21_fabs(INFINITY));
}
START_TEST(test_fabs_9) {
    ck_assert_double_infinite(s21_fabs(-INFINITY));
}
START_TEST(test_fabs_10) {
    ck_assert_double_nan(s21_fabs(NAN));
}
START_TEST(test_fabs_11) {
    double num = -3;
    ck_assert_double_eq_tol(s21_fabs(num), fabs(num), 1e-6);
}

START_TEST(test_fmax_1) {
    double num = 3;
    double y = 3;
    ck_assert_double_eq_tol(s21_fmax(num, y), fmax(num, y), 1e-6);
}
START_TEST(test_fmax_2) {
    double num = 3;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_fmax(num, y), fmax(num, y), 1e-6);
}
START_TEST(test_fmax_3) {
    double num = 100;
    double y = 5;
    ck_assert_double_eq_tol(s21_fmax(num, y), fmax(num, y), 1e-6);
}
START_TEST(test_fmax_4) {
    double num = 124;
    double y = 97;
    ck_assert_double_eq_tol(s21_fmax(num, y), fmax(num, y), 1e-6);
}
START_TEST(test_fmax_5) {
    double num = 3.1435;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_fmax(num, y), fmax(num, y), 1e-6);
}
START_TEST(test_fmax_6) {
    double num = 3.1435;
    double y = 8.2;
    ck_assert_double_eq_tol(s21_fmax(num, y), fmax(num, y), 1e-6);
}
START_TEST(test_fmax_7) {
    double num = 0;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_fmax(num, y), fmax(num, y), 1e-6);
}
START_TEST(test_fmax_8) {
    double num = -0;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_fmax(num, y), fmax(num, y), 1e-6);
}
START_TEST(test_fmax_9) {
    double num = 2148;
    double y = 1.2;
    ck_assert_double_eq_tol(s21_fmax(num, y), fmax(num, y), 1e-6);
}
START_TEST(test_fmax_10) {
    double num = NAN;
    double y = NAN;
    ck_assert_double_nan(s21_fmax(num, y));
}
START_TEST(test_fmax_12) {
    double num = NAN;
    double y = 1;
    ck_assert_double_eq_tol(s21_fmax(num, y), fmax(num, y), 1e-6);
    num = 1;
    y = NAN;
    ck_assert_double_eq_tol(s21_fmax(num, y), fmax(num, y), 1e-6);
}

Suite *s21_math_suite(void) {
    Suite *s;
    TCase *tc_core;
  TCase *tc_ceil;
  TCase *tc_floor;
  TCase *tc_tan;
  TCase *tc_cos;
  TCase *tc_sin;
  TCase *tc_abs;
  TCase *tc_atan;
  TCase *tc_acos;
  TCase *tc_asin;
  TCase *tc_fmod;
  TCase *tc_exp;
  TCase *tc_log;
  TCase *tc_sqrt;
  TCase *tc_pow;
  TCase *tc_fabs;
    s = suite_create("S21_math");
    tc_core = tcase_create("Core");

    // tcase_add_test(tc_core, test_abs);
    // tcase_add_test(tc_core, test_sin);
    // tcase_add_test(tc_core, test_sqrt);

    tcase_add_test(tc_core, test_log_1);
    tcase_add_test(tc_core, test_log_2);
    tcase_add_test(tc_core, test_log_3);
    tcase_add_test(tc_core, test_log_4);
    tcase_add_test(tc_core, test_log_5);
    tcase_add_test(tc_core, test_log_6);
    tcase_add_test(tc_core, test_log_7);
    tcase_add_test(tc_core, test_log_8);
    tcase_add_test(tc_core, test_log_9);
    tcase_add_test(tc_core, test_log_10);
    tcase_add_test(tc_core, test_log_11);
    tcase_add_test(tc_core, test_log_12);

    tcase_add_test(tc_core, test_pow_1);
    tcase_add_test(tc_core, test_pow_2);
    tcase_add_test(tc_core, test_pow_3);
    tcase_add_test(tc_core, test_pow_4);
    tcase_add_test(tc_core, test_pow_5);
    tcase_add_test(tc_core, test_pow_6);
    tcase_add_test(tc_core, test_pow_7);
    tcase_add_test(tc_core, test_pow_8);
    tcase_add_test(tc_core, test_pow_9);
    tcase_add_test(tc_core, test_pow_10);
    tcase_add_test(tc_core, test_pow_12);
    tcase_add_test(tc_core, test_pow_17);
    tcase_add_test(tc_core, test_pow_18);
    tcase_add_test(tc_core, test_pow_19);
    tcase_add_test(tc_core, test_pow_20);
    tcase_add_test(tc_core, test_pow_21);

    tcase_add_test(tc_core, test_fmod_1);
    tcase_add_test(tc_core, test_fmod_2);
    tcase_add_test(tc_core, test_fmod_3);
    tcase_add_test(tc_core, test_fmod_4);
    tcase_add_test(tc_core, test_fmod_5);
    tcase_add_test(tc_core, test_fmod_6);
    tcase_add_test(tc_core, test_fmod_7);
    tcase_add_test(tc_core, test_fmod_8);
    tcase_add_test(tc_core, test_fmod_9);
    tcase_add_test(tc_core, test_fmod_10);
    tcase_add_test(tc_core, test_fmod_11);
    tcase_add_test(tc_core, test_fmod_12);
    tcase_add_test(tc_core, test_fmod_14);
    tcase_add_test(tc_core, test_fmod_15);
    tcase_add_test(tc_core, test_fmod_16);
    tcase_add_test(tc_core, test_fmod_17);

    tcase_add_test(tc_core, test_fabs_1);
    tcase_add_test(tc_core, test_fabs_2);
    tcase_add_test(tc_core, test_fabs_3);
    tcase_add_test(tc_core, test_fabs_4);
    tcase_add_test(tc_core, test_fabs_5);
    tcase_add_test(tc_core, test_fabs_6);
    tcase_add_test(tc_core, test_fabs_7);
    tcase_add_test(tc_core, test_fabs_8);
    tcase_add_test(tc_core, test_fabs_9);
    tcase_add_test(tc_core, test_fabs_10);
    tcase_add_test(tc_core, test_fabs_11);

    tcase_add_test(tc_core, test_acos_1);
    tcase_add_test(tc_core, test_acos_2);
    tcase_add_test(tc_core, test_acos_3);
    tcase_add_test(tc_core, test_acos_4);
    tcase_add_test(tc_core, test_acos_5);
    tcase_add_test(tc_core, test_acos_6);
    tcase_add_test(tc_core, test_acos_7);
    tcase_add_test(tc_core, test_acos_8);
    tcase_add_test(tc_core, test_acos_9);
    tcase_add_test(tc_core, test_acos_10);
    tcase_add_test(tc_core, test_acos_11);
    tcase_add_test(tc_core, test_acos_12);
    tcase_add_test(tc_core, test_acos_13);
    tcase_add_test(tc_core, test_acos_14);
    tcase_add_test(tc_core, test_acos_15);

    tcase_add_test(tc_core, test_asin_1);
    tcase_add_test(tc_core, test_asin_2);
    tcase_add_test(tc_core, test_asin_3);
    tcase_add_test(tc_core, test_asin_4);
    tcase_add_test(tc_core, test_asin_5);
    tcase_add_test(tc_core, test_asin_6);
    tcase_add_test(tc_core, test_asin_7);
    tcase_add_test(tc_core, test_asin_8);
    tcase_add_test(tc_core, test_asin_9);
    tcase_add_test(tc_core, test_asin_10);
    tcase_add_test(tc_core, test_asin_11);

    tcase_add_test(tc_core, test_tan_1);
    tcase_add_test(tc_core, test_tan_2);
    tcase_add_test(tc_core, test_tan_3);
    tcase_add_test(tc_core, test_tan_4);
    tcase_add_test(tc_core, test_tan_5);
    tcase_add_test(tc_core, test_tan_6);
    tcase_add_test(tc_core, test_tan_7);
    tcase_add_test(tc_core, test_tan_8);
    tcase_add_test(tc_core, test_tan_9);
    tcase_add_test(tc_core, test_tan_10);
    tcase_add_test(tc_core, test_tan_11);
    tcase_add_test(tc_core, test_tan_12);

    tcase_add_test(tc_core, test_atan_1);
    tcase_add_test(tc_core, test_atan_2);
    tcase_add_test(tc_core, test_atan_3);
    tcase_add_test(tc_core, test_atan_4);
    tcase_add_test(tc_core, test_atan_5);
    tcase_add_test(tc_core, test_atan_6);
    tcase_add_test(tc_core, test_atan_7);
    tcase_add_test(tc_core, test_atan_8);
    tcase_add_test(tc_core, test_atan_9);
    tcase_add_test(tc_core, test_atan_10);
    tcase_add_test(tc_core, test_atan_11);

    tcase_add_test(tc_core, test_ceil_1);
    tcase_add_test(tc_core, test_ceil_2);
    tcase_add_test(tc_core, test_ceil_3);
    tcase_add_test(tc_core, test_ceil_4);
    tcase_add_test(tc_core, test_ceil_5);
    tcase_add_test(tc_core, test_ceil_6);
    tcase_add_test(tc_core, test_ceil_7);
    tcase_add_test(tc_core, test_ceil_8);
    tcase_add_test(tc_core, test_ceil_9);
    tcase_add_test(tc_core, test_ceil_10);
    tcase_add_test(tc_core, test_ceil_11);
    tcase_add_test(tc_core, test_ceil_12);
    tcase_add_test(tc_core, test_ceil_13);
    tcase_add_test(tc_core, test_ceil_14);
    tcase_add_test(tc_core, test_ceil_15);

    tcase_add_test(tc_core, test_floor_1);
    tcase_add_test(tc_core, test_floor_2);
    tcase_add_test(tc_core, test_floor_3);
    tcase_add_test(tc_core, test_floor_4);
    tcase_add_test(tc_core, test_floor_5);
    tcase_add_test(tc_core, test_floor_6);
    tcase_add_test(tc_core, test_floor_7);
    tcase_add_test(tc_core, test_floor_8);
    tcase_add_test(tc_core, test_floor_9);
    tcase_add_test(tc_core, test_floor_10);
    tcase_add_test(tc_core, test_floor_11);
    tcase_add_test(tc_core, test_floor_12);
    tcase_add_test(tc_core, test_floor_13);
    tcase_add_test(tc_core, test_floor_14);
    tcase_add_test(tc_core, test_floor_15);

    tcase_add_test(tc_core, test_cos_3);
    tcase_add_test(tc_core, test_cos_4);
    tcase_add_test(tc_core, test_cos_9);
    tcase_add_test(tc_core, test_cos_10);
    tcase_add_test(tc_core, test_cos_11);

    tcase_add_test(tc_core, test_exp_1);
    tcase_add_test(tc_core, test_exp_2);
    tcase_add_test(tc_core, test_exp_3);
    tcase_add_test(tc_core, test_exp_4);
    tcase_add_test(tc_core, test_exp_5);
    tcase_add_test(tc_core, test_exp_6);
    tcase_add_test(tc_core, test_exp_7);
    tcase_add_test(tc_core, test_exp_8);
    tcase_add_test(tc_core, test_exp_9);
    tcase_add_test(tc_core, test_exp_10);
    tcase_add_test(tc_core, test_exp_11);
    tcase_add_test(tc_core, test_exp_12);
    tcase_add_test(tc_core, test_exp_13);
    tcase_add_test(tc_core, test_exp_14);
    tcase_add_test(tc_core, test_exp_15);
    tcase_add_test(tc_core, test_exp_16);
    tcase_add_test(tc_core, test_exp_17);
    tcase_add_test(tc_core, test_exp_18);
    tcase_add_test(tc_core, test_exp_19);

    tcase_add_test(tc_core, test_fmax_1);
    tcase_add_test(tc_core, test_fmax_2);
    tcase_add_test(tc_core, test_fmax_3);
    tcase_add_test(tc_core, test_fmax_4);
    tcase_add_test(tc_core, test_fmax_5);
    tcase_add_test(tc_core, test_fmax_6);
    tcase_add_test(tc_core, test_fmax_7);
    tcase_add_test(tc_core, test_fmax_8);
    tcase_add_test(tc_core, test_fmax_9);
    tcase_add_test(tc_core, test_fmax_10);
    tcase_add_test(tc_core, test_fmax_12);
  tc_fabs = tcase_create("Fabs func");
  suite_add_tcase(s, tc_fabs);
  tcase_add_test(tc_fabs, s21_fabs_f);

  tc_pow = tcase_create("Pow func");
  suite_add_tcase(s, tc_pow);
  tcase_add_test(tc_pow, s21_pow_f);

  tc_sqrt = tcase_create("Sqrt func");
  suite_add_tcase(s, tc_sqrt);
  tcase_add_test(tc_sqrt, s21_sqrt_f);

  tc_log = tcase_create("Log func");
  suite_add_tcase(s, tc_log);
  tcase_add_test(tc_log, s21_log_f);

  tc_exp = tcase_create("Exp func");
  suite_add_tcase(s, tc_exp);
  tcase_add_test(tc_exp, s21_exp_f);

  tc_fmod = tcase_create("Fmod func");
  suite_add_tcase(s, tc_fmod);
  tcase_add_test(tc_fmod, s21_fmod_f);

  tc_asin = tcase_create("Asin func");
  suite_add_tcase(s, tc_asin);
  tcase_add_test(tc_asin, s21_asin_f);

  tc_acos = tcase_create("Acos func");
  suite_add_tcase(s, tc_acos);
  tcase_add_test(tc_acos, s21_acos_f);

  tc_atan = tcase_create("Atan func");
  suite_add_tcase(s, tc_atan);
  tcase_add_test(tc_atan, s21_atan_f);

  tc_abs = tcase_create("Abs func");
  suite_add_tcase(s, tc_abs);
  tcase_add_test(tc_abs, s21_abs_f);

  tc_sin = tcase_create("Sin func");
  suite_add_tcase(s, tc_sin);
  tcase_add_test(tc_sin, s21_sin_f);

  tc_cos = tcase_create("Cos func");
  suite_add_tcase(s, tc_cos);
  tcase_add_test(tc_cos, s21_cos_f);
tc_tan = tcase_create("Tan func");
  suite_add_tcase(s, tc_tan);
  tcase_add_test(tc_tan, s21_tan_f);

  tc_floor = tcase_create("Floor func");
  suite_add_tcase(s, tc_floor);
  tcase_add_test(tc_floor, s21_floor_f);

  tc_ceil = tcase_create("Ceil func");
  suite_add_tcase(s, tc_ceil);
  tcase_add_test(tc_ceil, s21_ceil_f);

    suite_add_tcase(s, tc_core);
    return s;
}

int main(void) {
    Suite *s;
    SRunner *sr;

    s = s21_math_suite();
    sr = srunner_create(s);

    srunner_set_fork_status(sr, CK_NOFORK);
    srunner_run_all(sr, CK_NORMAL);
    srunner_free(sr);
    return 0;
}
