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

#endif // DSP_OPTIMIZATION
