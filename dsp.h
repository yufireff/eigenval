#ifndef _DSP_H_
#define _DSP_H_

extern int RealMatrixMult;
extern int ComplexMatrixMult;
extern int ComplexMatrixInv;
extern int RealMatrixMultLeftTransp;
extern int ComplexMatrixMultLeftTransp;
extern int ComplexMatrixMultRightTransp;
extern int RealMatrixMultAtBA;
extern int MatrixAminusUSUt;
extern int QRSymmReal;

extern float g_Buffer1[49*49];
extern float g_Buffer2[49*49];
extern float g_Buffer3d[49*49];
extern float g_Buffer4d[49*49];
extern float g_Buffer5[49*49];
extern float g_Buffer6[49*49];
extern float g_Buffer7[49*49];
extern float g_pDiag[49];

extern float g_Factor;
extern float g_FactorIm;
extern int g_nRows1;
extern int g_nColumns1;
extern int g_nRows2;
extern int g_nColumns2;
extern int g_Size;
extern float g_Norm;

void init_dsp(int dspNum);
void run_dsp(int dspNum, unsigned int command);

#define DSP_ROUTINE_ADDR(x) ((unsigned int)&x - 0xb8440000)>>2

unsigned int GetCP0_Count();

#endif // _DSP_H_
