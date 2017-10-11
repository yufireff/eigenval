#include "matrix.h"
#include "MatrixSpecial.h"

int complex_matrix_mult_right_transp(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res)
{
	int i, j, k, i1, j1;

	REAL_TYPE real, img;
	REAL_TYPE x, y;

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
				i1 = k*a->numCols + j;
				j1 = k*b->numCols + j;
				sum += a->pData[i1] * b->pData[j1];
			}
			res->pData[i*res->numCols + j] = factor * sum;
		}
	}

	return MATRIX_SUCCESS;
}