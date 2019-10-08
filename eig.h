#ifndef EIG_H
#define EIG_H

#include "matrix.h"

void eig(const CMatrix_t* A, CMatrix_t* e, CMatrix_t* ev);

void complex_triag(const CMatrix_t* A, Matrix_t* T, CMatrix_t* P);

void qr_symm_real(const Matrix_t* T, REAL_TYPE tol, Matrix_t* U, Matrix_t* D);

int eig_symm_triag(const CMatrix_t* A, REAL_TYPE tol, Matrix_t* S, CMatrix_t* U);

REAL_TYPE eig_symm_triag_only_one(const CMatrix_t* A, REAL_TYPE tol, CMatrix_t* U);

#endif // EIG_H
