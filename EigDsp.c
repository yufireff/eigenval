#include "EigDsp.h"
#include "eig.h"
#include "dsp.h"

#include <string.h>

#ifdef DSP_OPTIMIZATION_FULL
void qr_symm_real_dsp(const Matrix_t* T, REAL_TYPE tol, Matrix_t* U, Matrix_t* D, int numDsp)
{
	init_dsp(numDsp);

        memcpy(&(g_Buffer1[0]), T->pData, T->numRows * T->numCols * sizeof(float));
        memcpy(&g_Factor, &tol, sizeof(float));
        memcpy(&g_Size, &U->numRows, sizeof(int));

        run_dsp(numDsp, DSP_ROUTINE_ADDR(QRSymmReal));
	memcpy(U->pData, &(g_Buffer7[0]), U->numRows * U->numCols * sizeof(float));
	memcpy(D->pData, &(g_Buffer1[0]), D->numRows * D->numCols * sizeof(float));
}
#else
void qr_symm_real_dsp(const Matrix_t* T, REAL_TYPE tol, Matrix_t* U, Matrix_t* D, int numDsp)
{
	qr_symm_real(T, tol, U, D);
}
#endif // DSP_OPTIMIZATION_FULL
