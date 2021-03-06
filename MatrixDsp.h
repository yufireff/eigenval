#ifndef _MATRIX_DSP_H_
#define _MATRIX_DSP_H_

#include "settings.h"

#ifdef DSP_OPTIMIZATION

#ifdef DOUBLE
#pragma error "DSP functions only for single float"
#endif // DOUBLE

#endif // DSP_OPTIMIZATION

#include "matrix.h"

#ifndef DSP_OPTIMIZATION_FULL
int real_matrix_mult_dsp(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor_re, Matrix_t* res, int numDsp);

int complex_matrix_mult_dsp(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res, int numDsp);

int inv_dsp(const CMatrix_t* a, CMatrix_t* ainv, int numDsp);

int real_matrix_mult_left_transp_dsp(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor_re, Matrix_t* res, int numDsp);
#endif // !DSP_OPTIMIZATION_FULL

int complex_matrix_mult_left_transp_dsp(const CMatrix_t* a, const CMatrix_t* b, CMatrix_t* res, int numDsp);

#ifndef DSP_OPTIMIZATION_FULL
int complex_matrix_mult_right_transp_dsp(const CMatrix_t* a, const CMatrix_t* b, CMatrix_t* res, int numDsp);

int real_matrix_mult_at_b_a_dsp(const Matrix_t* a, const Matrix_t* b, Matrix_t* res, int numDsp);
#endif // !DSP_OPTIMIZATION_FULL

#ifdef TEST_RESULT
int complex_a_minus_u_s_ut_dsp(const CMatrix_t* a, const CMatrix_t* u, const Matrix_t* s, CMatrix_t *res, float *norm, int numDsp);
#endif // TEST_RESULT
#endif // _MATRIX_DSP_H_
