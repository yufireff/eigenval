 .global g_pReal1
 .global g_pImag1
 .global g_pReal2
 .global g_pImag2
 .global g_pReal3
 .global g_pImag3
 .global g_Factor
 .global g_FactorIm
 .global g_nRows1
 .global g_nColumns1
 .global g_nColumns2

 .data
 g_pReal1: .space 9604, 0
 g_pImag1: .space 9604, 0
 g_pReal2: .space 19208, 0
 g_pImag2: .space 19208, 0
 g_pReal3: .space 9604, 0
 g_pImag3: .space 9604, 0
 g_Factor: .real 0
 g_nRows1: .word 0
 g_nColumns1: .word 0
 g_nColumns2: .word 0
 g_FactorIm: .real 0