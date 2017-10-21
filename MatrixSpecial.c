#include "matrix.h"
#include "MatrixSpecial.h"

#include "complex.h"

#ifdef PREALLOCATION
#ifdef TEST_RESULT
static REAL_TYPE buffer[MATRIX_MAX_SIZE * 6];
#else // TEST_RESULT
static REAL_TYPE buffer[MATRIX_MAX_SIZE];
#endif // TEST_RESULT
#endif // PREALLOCATION

int complex_matrix_mult_right_transp(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res)
{
	int i, j, k, i1, j1;

	REAL_TYPE real, img;
	REAL_TYPE x = 0.0f, y = 0.0f;

	if (a->numCols != b->numCols || res->numRows != a->numRows || res->numCols != b->numRows)
		return MATRIX_EXCEEDS;

	for (i = 0; i < a->numRows; ++i)
	{
		for (j = 0; j < b->numRows; ++j)
		{
			real = 0.0;
			img = 0.0;
			for (k = 0; k < a->numCols; ++k)
			{
				i1 = i*a->numCols + k;
				j1 = j*b->numCols + k;
				complex_mult(a->pDataReal[i1], a->pDataImag[i1],
					b->pDataReal[j1], -b->pDataImag[j1], &x, &y);
				real += x;
				img += y;
			}
			k = i*res->numCols + j;
			complex_mult(real, img, factor_re, factor_im, res->pDataReal + k, res->pDataImag + k);
		}
	}

	return MATRIX_SUCCESS;
}

int real_matrix_mult_left_transp(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor, Matrix_t* res)
{
	int i, j, k, i1, j1;

	REAL_TYPE sum;

	if (a->numRows != b->numRows || res->numRows != a->numCols || res->numCols != b->numCols)
		return MATRIX_EXCEEDS;

	for (i  = 0; i < a->numCols; ++i)
	{
		for (j = 0; j < b->numCols; ++j)
		{
			sum = 0.0f;
			for (k = 0; k < a->numRows; ++k)
			{
				i1 = k*a->numCols + i;
				j1 = k*b->numCols + j;
				sum += a->pData[i1] * b->pData[j1];
			}
			res->pData[i*res->numCols + j] = factor * sum;
		}
	}

	return MATRIX_SUCCESS;
}

int complex_matrix_mult_right_transp_d(const CMatrix_td* a, const CMatrix_td* b, CMatrix_td* res)
{
	int i, j, k;

	double real, imag, x = 0.0, y = 0.0;
	double *pReal1, *pImag1;
	double *pReal2, *pImag2;
	double *pReal3 = res->pDataReal, *pImag3 = res->pDataImag;
	
	if (a->numCols != b->numCols || res->numRows != a->numRows || res->numCols != b->numRows)
		return MATRIX_EXCEEDS;

	for (i = 0; i < a->numRows; ++i)
	{
		for (j = 0; j < b->numRows; ++j)
		{
			k = i * a->numCols;
			pReal1 = a->pDataReal + k;
			pImag1 = a->pDataImag + k;
			k = j * b->numCols;
			pReal2 = b->pDataReal + k;
			pImag2 = b->pDataImag + k;
			real = 0.0;
			imag = 0.0;
			for (k = 0; k < a->numCols; ++k)
			{
				complex_mult_d(*pReal1++, *pImag1++,
					*pReal2++, -(*pImag2++), &x, &y);
				real += x;
				imag += y;
			}
			*pReal3++ = real;
			*pImag3++ = imag;
		}
	}

	return MATRIX_SUCCESS;
}

int complex_matrix_mult_left_transp_d(const CMatrix_td* a, const CMatrix_td* b, CMatrix_td* res)
{
	int i, j, k, i1, j1;

	double real, imag;
	double x = 0.0, y = 0.0;

	if (a->numCols != b->numCols || res->numRows != a->numRows || res->numCols != b->numRows)
		return MATRIX_EXCEEDS;

	for (i = 0; i < a->numCols; ++i)
	{
		for (j = 0; j < b->numCols; ++j)
		{
			real = 0.0;
			imag = 0.0;
			for (k = 0; k < a->numCols; ++k)
			{
				i1 = k*a->numCols + i;
				j1 = k*b->numCols + j;
				complex_mult_d(a->pDataReal[i1], a->pDataImag[i1],
					-b->pDataReal[j1], b->pDataImag[j1], &x, &y);
				real += x;
				imag += y;
			}
			k = i*res->numCols + j;
			res->pDataReal[k] = real;
			res->pDataImag[k] = imag;
			//complex_mult_d(real, img, factor_re, factor_im, res->pDataReal + k, res->pDataImag + k);
		}
	}

	return MATRIX_SUCCESS;
}

int complex_scal_prod_left_transp_d(const CMatrix_td* a, const CMatrix_td* b, double* re, double* im)
{
	int i;
	double tmp_re = 0.0, tmp_im = 0.0;

	/*if (a->numRows != 1 || b->numCols != 1 || a->numCols != b->numRows)
	return MATRIX_EXCEEDS;*/

	*re = 0.0f;
	*im = 0.0f;

	for (i = 0; i < int_max(a->numRows, a->numCols); ++i)
	{
		complex_mult_d(a->pDataReal[i], -a->pDataImag[i], b->pDataReal[i], b->pDataImag[i], &tmp_re, &tmp_im);
		*re += tmp_re;
		*im += tmp_im;
	}

	return MATRIX_SUCCESS;
}

int complex_matrix_mult_left_transp(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res)
{
	int i, j, k, i1, j1;

	REAL_TYPE real, img;
	REAL_TYPE x = 0.0f, y = 0.0f;
	if (a->numRows != b->numRows || res->numRows != a->numCols || res->numCols != b->numCols)
		return MATRIX_EXCEEDS;

	for (i = 0; i < a->numCols; ++i)
	{
		for (j = 0; j < b->numCols; ++j)
		{
			real = 0.0;
			img = 0.0;
			for (k = 0; k < a->numRows; ++k)
			{
				i1 = k*a->numCols + i;
				j1 = k*b->numCols + j;
				complex_mult(a->pDataReal[i1], a->pDataImag[i1],
					b->pDataReal[j1], b->pDataImag[j1], &x, &y);
				real += x;
				img += y;
			}
			k = i*res->numCols + j;
			complex_mult(real, img, factor_re, factor_im, res->pDataReal + k, res->pDataImag + k);
		}
	}

	return MATRIX_SUCCESS;
}

int real_matrix_mult_at_b_a(const Matrix_t* a, const Matrix_t* b, Matrix_t* res)
{
	Matrix_t ba;
#ifdef PREALLOCATION
	ba.pData = buffer;
#endif // PREALLOCATION

	if (a->numRows != a->numCols || b->numRows != b->numCols || res->numCols != res->numRows
		|| a->numRows != b->numRows || a->numRows != res->numRows)
		return MATRIX_EXCEEDS;

	real_new(a->numRows, a->numCols, &ba);

	real_matrix_mult(b, a, 1.0f, &ba);
	real_matrix_mult_left_transp(a, &ba, 1.0f, res);

	return MATRIX_SUCCESS;
}

#ifdef TEST_RESULT
int complex_a_minus_u_s_ut(const CMatrix_t* a, const CMatrix_t* u, const Matrix_t* s, CMatrix_t *res, float *norm)
{
	CMatrix_t us, usut, scmpl;
	int n = a->numRows, i;
	if (a->numRows != a->numCols || u->numRows != u->numCols || res->numRows != res->numCols
		|| a->numRows != u->numRows || a->numRows != res->numRows || s->numRows != 1 || s->numCols != n)
		return MATRIX_EXCEEDS;

#ifdef PREALLOCATION
	us.pDataReal = buffer;
	us.pDataImag = buffer + MATRIX_MAX_SIZE;
	usut.pDataReal = buffer + 2 * MATRIX_MAX_SIZE;
	usut.pDataImag = buffer + 3 * MATRIX_MAX_SIZE;
	scmpl.pDataReal = buffer + 4 * MATRIX_MAX_SIZE;
	scmpl.pDataImag = buffer + 5 * MATRIX_MAX_SIZE;
#endif // PREALLOCATION
	complex_new(n, n, &us);
	complex_new(n, n, &usut);
	complex_zeros(n, n, &scmpl, 1);

	for (i = 0; i < n; ++i)
		scmpl.pDataReal[i*n + i] = s->pData[i];

	complex_matrix_mult(u, &scmpl, 1.0f, 0.0f, &us);
	complex_matrix_mult_right_transp(&us, u, 1.0f, 0.0f, &usut);
	complex_matrix_sum(&usut, a, -1.0f, 0.0f, res);

	*norm = 0.0f;
	for (i = 0; i < n*n; ++i)
	{
		*norm += res->pDataReal[i] * res->pDataReal[i];
		*norm += res->pDataImag[i] * res->pDataImag[i];
	}

	complex_free(&us);
	complex_free(&usut);
	complex_free(&scmpl);
	return MATRIX_SUCCESS;
}
#endif // TEST_RESULT

