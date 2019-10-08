#ifndef EIG_DSP_H
#define EIG_DSP_H

#include "matrix.h"

void qr_symm_real_dsp(const Matrix_t* T, REAL_TYPE tol, Matrix_t* U, Matrix_t* D, int numDsp);

#endif // EIG_DSP_H
