#include "MatrixDsp.h"

#ifdef DSP_OPTIMIZATION

#include <string.h>
#include "dsp.h"

int complex_matrix_mult_dsp(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res, int numDsp)
{
	init_dsp(numDsp);

	memcpy(&(g_pReal1[0]), a->pDataReal, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pReal2[0]), b->pDataReal, b->numRows * b->numCols * sizeof(float));
	memcpy(&(g_pImag1[0]), a->pDataImag, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pImag2[0]), b->pDataImag, b->numRows * b->numCols * sizeof(float));

	g_Factor = factor_re;
	g_FactorIm = factor_im;
	g_nRows1 = a->numRows;
	g_nColumns1 = a->numCols;
	g_nColumns2 = b->numCols;

	run_dsp(numDsp, DSP_ROUTINE_ADDR(ComplexMatrixMult));
	memcpy(res->pDataReal, &(g_pReal3[0]), res->numRows * res->numCols * sizeof(float));
	memcpy(res->pDataImag, &(g_pImag3[0]), res->numRows * res->numCols * sizeof(float));

	return MATRIX_SUCCESS;
}

int real_matrix_mult_dsp(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor, Matrix_t* res, int numDsp)
{
	init_dsp(numDsp);

	memcpy(&(g_pReal1[0]), a->pData, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pReal2[0]), b->pData, b->numRows * b->numCols * sizeof(float));

	g_Factor = factor;
	g_nRows1 = a->numRows;
	g_nColumns1 = a->numCols;
	g_nColumns2 = b->numCols;

	run_dsp(numDsp, DSP_ROUTINE_ADDR(RealMatrixMult));
	memcpy(res->pData, &(g_pReal3[0]), res->numRows * res->numCols * sizeof(float));

	return MATRIX_SUCCESS;
}

int inv_dsp(const CMatrix_t* a, CMatrix_t* ainv, int numDsp)
{
    init_dsp(numDsp);

	memcpy(&(g_pReal1[0]), a->pDataReal, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pImag1[0]), a->pDataImag, a->numRows * a->numCols * sizeof(float));

	g_nRows1 = a->numRows;
	run_dsp(numDsp, DSP_ROUTINE_ADDR(ComplexMatrixInv));

	memcpy(ainv->pDataReal, &(g_pReal1[0]), ainv->numRows * ainv->numCols * sizeof(float));
	memcpy(ainv->pDataImag, &(g_pImag1[0]), ainv->numRows * ainv->numCols * sizeof(float));
	return MATRIX_SUCCESS;

}

int real_matrix_mult_left_transp_dsp(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor, Matrix_t* res, int numDsp)
{
	init_dsp(numDsp);

	memcpy(&(g_pReal1[0]), a->pData, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pReal2[0]), b->pData, b->numRows * b->numCols * sizeof(float));

	g_Factor = factor;
	g_nRows1 = a->numRows;
	g_nColumns1 = a->numCols;
	g_nColumns2 = b->numCols;

	run_dsp(numDsp, DSP_ROUTINE_ADDR(RealMatrixMultLeftTransp));
	memcpy(res->pData, &(g_pReal3[0]), res->numRows * res->numCols * sizeof(float));

	return MATRIX_SUCCESS;
}

int complex_matrix_mult_left_transp_dsp(const CMatrix_t* a, const CMatrix_t* b, CMatrix_t* res, int numDsp)
{
	init_dsp(numDsp);

	memcpy(&(g_pReal1[0]), a->pDataReal, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pReal2[0]), b->pDataReal, b->numRows * b->numCols * sizeof(float));
	memcpy(&(g_pImag1[0]), a->pDataImag, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pImag2[0]), b->pDataImag, b->numRows * b->numCols * sizeof(float));


	g_nRows1 = a->numRows;
	g_nColumns1 = a->numCols;
	g_nColumns2 = b->numCols;

	run_dsp(numDsp, DSP_ROUTINE_ADDR(ComplexMatrixMultLeftTransp));
	memcpy(res->pDataReal, &(g_pReal3[0]), res->numRows * res->numCols * sizeof(float));
	memcpy(res->pDataImag, &(g_pImag3[0]), res->numRows * res->numCols * sizeof(float));

	return MATRIX_SUCCESS;
}

int complex_matrix_mult_right_transp_dsp(const CMatrix_t* a, const CMatrix_t* b, CMatrix_t* res, int numDsp)
{
	init_dsp(numDsp);

	memcpy(&(g_pReal1[0]), a->pDataReal, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pReal2[0]), b->pDataReal, b->numRows * b->numCols * sizeof(float));
	memcpy(&(g_pImag1[0]), a->pDataImag, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pImag2[0]), b->pDataImag, b->numRows * b->numCols * sizeof(float));


	g_nRows1 = a->numRows;
	g_nColumns1 = a->numCols;
	g_nRows2 = b->numRows;

	run_dsp(numDsp, DSP_ROUTINE_ADDR(ComplexMatrixMultRightTransp));
	memcpy(res->pDataReal, &(g_pReal3[0]), res->numRows * res->numCols * sizeof(float));
	memcpy(res->pDataImag, &(g_pImag3[0]), res->numRows * res->numCols * sizeof(float));

	return MATRIX_SUCCESS;
}

int real_matrix_mult_at_b_a_dsp(const Matrix_t* a, const Matrix_t* b, Matrix_t* res, int numDsp)
{
	init_dsp(numDsp);

	memcpy(&(g_pReal1[0]), a->pData, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pReal2[0]), b->pData, b->numRows * b->numCols * sizeof(float));

	g_Size = a->numRows;

	run_dsp(numDsp, DSP_ROUTINE_ADDR(RealMatrixMultAtBA));
	memcpy(res->pData, &(g_pReal3[0]), res->numRows * res->numCols * sizeof(float));

	return MATRIX_SUCCESS;
}

int complex_a_minus_u_s_ut_dsp(const CMatrix_t* a, const CMatrix_t* u, const Matrix_t* s, CMatrix_t *res, float *norm, int numDsp)
{
	init_dsp(numDsp);

	memcpy(&(g_pReal1[0]), u->pDataReal, u->numRows * u->numCols * sizeof(float));
	memcpy(&(g_pImag1[0]), u->pDataImag, u->numRows * u->numCols * sizeof(float));
	memcpy(&(g_pReal2[0]), a->pDataReal, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pImag2[0]), a->pDataImag, a->numRows * a->numCols * sizeof(float));
	memcpy(&(g_pDiag[0]), s->pData, s->numRows * s->numCols * sizeof(float));

	g_Size = a->numRows;

	run_dsp(numDsp, DSP_ROUTINE_ADDR(MatrixAminusUSUt));

	memcpy(res->pDataReal, &(g_pReal3[0]), res->numRows * res->numCols * sizeof(float));
	memcpy(res->pDataImag, &(g_pImag3[0]), res->numRows * res->numCols * sizeof(float));
	*norm = g_Norm;

	return MATRIX_SUCCESS;
}

#else // DSP_OPTIMIZATION

#include "matrix.h"
#include "MatrixSpecial.h"

int complex_matrix_mult_dsp(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res, int numDsp)
{
	CMatrix_t copy;
#ifdef PREALLOCATION
	static REAL_TYPE buffer[2 * MATRIX_MAX_SIZE];
	copy.pDataReal = buffer;
	copy.pDataImag = buffer + MATRIX_MAX_SIZE;
#endif // PREALLOCATION
	if (a == res)
	{
		complex_clone(a, &copy, 1);
		return complex_matrix_mult(&copy, b, factor_re, factor_im, res);
	}
	else if (b == res)
	{
		complex_clone(b, &copy, 1);
		return complex_matrix_mult(a, &copy, factor_re, factor_im, res);
	}
	else
		return complex_matrix_mult(a, b, factor_re, factor_im, res);
}

int real_matrix_mult_dsp(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor, Matrix_t* res, int numDsp)
{
	return real_matrix_mult(a, b, factor, res);
}

int inv_dsp(const CMatrix_t* a, CMatrix_t* ainv, int numDsp)
{
	return inv(a, ainv);
}

int real_matrix_mult_left_transp_dsp(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor, Matrix_t* res, int numDsp)
{
	return real_matrix_mult_left_transp(a, b, factor, res);
}

int complex_matrix_mult_left_transp_dsp(const CMatrix_t* a, const CMatrix_t* b, CMatrix_t* res, int numDsp)
{
	return complex_matrix_mult_left_transp(a, b, 1.0f, 0.0f, res);
}

int complex_matrix_mult_right_transp_dsp(const CMatrix_t* a, const CMatrix_t* b, CMatrix_t* res, int numDsp)
{
	return complex_matrix_mult_right_transp(a, b, 1.0f, 0.0f, res);
}

int real_matrix_mult_at_b_a_dsp(const Matrix_t* a, const Matrix_t* b, Matrix_t* res, int numDsp)
{
	return real_matrix_mult_at_b_a(a, b, res);
}

int complex_a_minus_u_s_ut_dsp(const CMatrix_t* a, const CMatrix_t* u, const Matrix_t* s, CMatrix_t *res, float *norm, int numDsp)
{
	return complex_a_minus_u_s_ut(a, u, s, res, norm);
}
#endif // DSP_OPTIMIZATION
