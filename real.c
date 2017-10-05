#include "real.h"
#include <float.h>
#include <math.h>

int int_max(int a, int b)
{
	return (a > b) ? a : b;
}

REAL_TYPE sign(REAL_TYPE x)
{
	return (x < 0.0f) ? -1.0f : 1.0f;
}

double sign_d(double x)
{
	return (x < 0.0) ? -1.0 : 1.0;
}

int is_zero(REAL_TYPE x)
{
#ifdef DOUBLE
	return (fabs(x) <= DBL_EPSILON) ? 1 : 0;
#else
	return (fabsà(x) <= FLT_EPSILON) ? 1 : 0;
#endif
	
}

int is_zero_d(double x)
{
	return (fabs(x) <= DBL_EPSILON) ? 1 : 0;
}

