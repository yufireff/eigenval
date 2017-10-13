#ifndef _DSP_H_
#define _DSP_H_

extern int RealMatrixMult;
extern int ComplexMatrixMult;
extern int ComplexMatrixInv;
extern int RealMatrixMultLeftTransp;
extern int ComplexMatrixMultLeftTransp;
extern int ComplexMatrixMultRightTransp;

extern float g_pReal1[49*49];
extern float g_pImag1[49*49];
extern float g_pReal2[49*49];
extern float g_pImag2[49*49];
extern float g_pReal3[49*49];
extern float g_pImag3[49*49];

extern float g_Factor;
extern float g_FactorIm;
extern int g_nRows1;
extern int g_nColumns1;
extern int g_nRows2;
extern int g_nColumns2;

void init_dsp(int dspNum);
void run_dsp(int dspNum, unsigned int command);

#define DSP_ROUTINE_ADDR(x) ((unsigned int)&x - 0xb8440000)>>2

#endif // _DSP_H_
