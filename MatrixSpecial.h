#ifndef _MATRIX_SPECIAL_H_
#define _MATRIX_SPECIAL_H_

int complex_matrix_mult_right_transp(const CMatrix_t* a, const CMatrix_t* b, REAL_TYPE factor_re, REAL_TYPE factor_im, CMatrix_t* res);

int real_matrix_mult_left_transp(const Matrix_t* a, const Matrix_t* b, REAL_TYPE factor, Matrix_t* res);


#endif //
