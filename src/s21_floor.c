#include "s21_math.h"

long double s21_floor(double x) {
  double arg = x;
  int num = arg;
  long double res;
  if (x == S21_INFINITY_POSITIVE || x == S21_INFINITY_NEGATIVE) {
    res = x;
  } else if (arg == num) {
    res = num;
  } else if (arg > (long double)num) {
    res = num;
  } else if (arg < (long double)num) {
    res = num - 1;
  }

  return x != x ? S21_NAN : res;
}
