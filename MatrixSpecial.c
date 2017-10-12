#include "matrix.h"
#include "MatrixSpecial.h"

#include "complex.h"

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
	int i, j, k, i1, j1;

	double real, imag;
	double x = 0.0, y = 0.0;

	if (a->numCols != b->numCols || res->numRows != a->numRows || res->numCols != b->numRows)
		return MATRIX_EXCEEDS;

	for (i = 0; i < a->numRows; ++i)
	{
		for (j = 0; j < b->numRows; ++j)
		{
			real = 0.0;
			imag = 0.0;
			for (k = 0; k < a->numCols; ++k)
			{
				i1 = i*a->numCols + k;
				j1 = j*b->numCols + k;
				complex_mult_d(a->pDataReal[i1], a->pDataImag[i1],
					b->pDataReal[j1], -b->pDataImag[j1], &x, &y);
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
