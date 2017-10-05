#ifndef COMPLEX_H
#define COMPLEX_H

#include "settings.h"

void complex_mult(REAL_TYPE r1, REAL_TYPE i1, REAL_TYPE r2, REAL_TYPE i2, REAL_TYPE* r, REAL_TYPE* i);

void complex_mult_d(double r1, double i1, double r2, double i2, double* r, double* i);

void complex_div(REAL_TYPE r1, REAL_TYPE i1, REAL_TYPE r2, REAL_TYPE i2, REAL_TYPE* r, REAL_TYPE* i);

void complex_div_d(double r1, double i1, double r2, double i2, double* r, double* i);

REAL_TYPE complex_sign(REAL_TYPE re, REAL_TYPE im);

double complex_sign_d(double re, double im);

int is_zero_complex(REAL_TYPE re, REAL_TYPE im);

REAL_TYPE complex_norm_square(REAL_TYPE re, REAL_TYPE im);

REAL_TYPE complex_norm(REAL_TYPE re, REAL_TYPE im);

#endif // COMPLEX_H
