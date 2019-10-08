#ifndef EIG_H
#define EIG_H

#include "matrix.h"

void complex_triag(const CMatrix_t* A, Matrix_t* T, CMatrix_t* P);
#ifndef DSP_OPTIMIZATION_FULL
void qr_symm_real(const Matrix_t* T, REAL_TYPE tol, Matrix_t* U, Matrix_t* D);
#endif // !DSP_OPTIMIZATION_FULL
int eig_symm_triag(const CMatrix_t* A, REAL_TYPE tol, Matrix_t* S, CMatrix_t* U);

#endif // EIG_H
