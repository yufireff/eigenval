#ifndef REAL_H
#define REAL_H

#include "settings.h"

int int_max(int a, int b);

REAL_TYPE sign(REAL_TYPE x);

double sign_d(double x);

int is_zero(REAL_TYPE x);

int is_zero_d(double x);

#ifdef DOUBLE
#define sq_rt sqrt
#define f_abs fabs
#else // DOUBLE
#define sq_rt sqrtf
#define f_abs fabsf
#endif

#endif // REAL_H
