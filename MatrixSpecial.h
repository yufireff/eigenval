#ifndef _MATRIX_SPECIAL_H_
#define _MATRIX_SPECIAL_H_

int complex_matrix_mult_right_transp(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res);

int real_matrix_mult_left_transp(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor, Matrix_t* res);

int complex_matrix_mult_left_transp_d(const CMatrix_td* a, const CMatrix_td* b,  CMatrix_td* res);

int complex_matrix_mult_right_transp_d(const CMatrix_td* a, const CMatrix_td* b, CMatrix_td* res);

int complex_scal_prod_left_transp_d(const CMatrix_td* a, const CMatrix_td* b, double* re, double* im);

int complex_matrix_mult_left_transp(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res);

int real_matrix_mult_at_b_a(const Matrix_t* a, const Matrix_t* b, Matrix_t* res);

int complex_a_minus_u_s_ut(const CMatrix_t* a, const CMatrix_t* u, const Matrix_t* s, CMatrix_t *res, float *norm);

#endif //
