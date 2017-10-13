#ifndef _MATRIX_DSP_H_
#define _MATRIX_DSP_H_

#include "settings.h"

#ifdef DSP_OPTIMIZATION

#ifdef DOUBLE
#pragma error "DSP functions only for single float"
#endif // DOUBLE

#endif // DSP_OPTIMIZATION

#include "matrix.h"

int real_matrix_mult_dsp(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor_re, Matrix_t* res, int numDsp);

int complex_matrix_mult_dsp(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res, int numDsp);

int inv_dsp(const CMatrix_t* a, CMatrix_t* ainv, int numDsp);

int real_matrix_mult_left_transp_dsp(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor_re, Matrix_t* res, int numDsp);

int complex_matrix_mult_left_transp_dsp(const CMatrix_t* a, const CMatrix_t* b, CMatrix_t* res, int numDsp);

int complex_matrix_mult_right_transp_dsp(const CMatrix_t* a, const CMatrix_t* b, CMatrix_t* res, int numDsp);

#endif // _MATRIX_DSP_H_
