#include "complex.h"

#include <math.h>

#include "real.h"


void complex_mult(REAL_TYPE r1, REAL_TYPE i1, REAL_TYPE r2, REAL_TYPE i2, REAL_TYPE* r, REAL_TYPE* i)
{
	*r = r1 * r2 - i1 * i2;
	*i = r1 * i2 + r2 * i1;
}

void complex_mult_d(double r1, double i1, double r2, double i2, double* r, double* i)
{
	*r = r1 * r2 - i1 * i2;
	*i = r1 * i2 + r2 * i1;
}

void complex_div(REAL_TYPE r1, REAL_TYPE i1, REAL_TYPE r2, REAL_TYPE i2, REAL_TYPE* r, REAL_TYPE* i)
{
	REAL_TYPE norm2 = r2 * r2 + i2 * i2;
	*r = (r1 * r2 + i1 * i2) / norm2;
	*i = (r2 * i1 - r1 * i2) / norm2;
}

void complex_div_d(double r1, double i1, double r2, double i2, double* r, double* i)
{
	double norm2 = r2 * r2 + i2 * i2;
	*r = (r1 * r2 + i1 * i2) / norm2;
	*i = (r2 * i1 - r1 * i2) / norm2;
}

REAL_TYPE complex_sign(REAL_TYPE re, REAL_TYPE im)
{
	if (!is_zero(re))
		return sign(re);
	else if (!is_zero(im))
		return sign(im);
	else
		return 1.0f;
}

double complex_sign_d(double re, double im)
{
	if (!is_zero_d(re))
		return sign_d(re);
	else if (!is_zero_d(im))
		return sign_d(im);
	else
		return 1.0;
}

int is_zero_complex(REAL_TYPE re, REAL_TYPE im)
{
	return (is_zero(re) && is_zero(im)) ? 0 : 1;
}

REAL_TYPE complex_norm_square(REAL_TYPE re, REAL_TYPE im)
{
	return re * re + im * im;
}

REAL_TYPE complex_norm(REAL_TYPE re, REAL_TYPE im)
{
	return sq_rt(re*re + im*im);
}

