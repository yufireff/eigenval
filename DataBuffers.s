 .global g_pReal1
 .global g_pImag1
 .global g_pReal2
 .global g_pImag2
 .global g_pReal3
 .global g_pImag3
 .global g_pReal4
 .global g_pImag4
 .global g_pDiag
 .global g_Factor
 .global g_FactorIm
 .global g_nRows1
 .global g_nColumns1
 .global g_nRows2
 .global g_nColumns2
 .global g_Size
 .global g_Norm

 .data
 g_pReal1: .space 49*49*4, 0
 g_pImag1: .space 49*49*4, 0
 g_pReal2: .space 2*49*49*4, 0
 g_pImag2: .space 2*49*49*4, 0
 g_pReal3: .space 49*49*4, 0
 g_pImag3: .space 49*49*4, 0
 g_pReal4: .space 49*49*4, 0
 g_pImag4: .space 49*49*4, 0
 g_pDiag:  .space 49*4, 0
 ; для отладки
 ;g_pReal1: .space 4*4*4, 0
 ;g_pImag1: .space 4*4*4, 0
 ;g_pReal2: .space 4*4*4, 0
 ;g_pImag2: .space 4*4*4, 0
 ;g_pReal3: .space 4*4*4, 0
 ;g_pImag3: .space 4*4*4, 0
 ;g_pDiag:  .space 4*4, 0
 ;g_pReal4: .space 4*4*4, 0
 ;g_pImag4: .space 4*4*4, 0

 g_Factor: .real 0
 g_nRows1: .word 0
 g_nColumns1: .word 0
 g_nRows2: .word 0
 g_nColumns2: .word 0
 g_FactorIm: .real 0
 g_Size: .word 0
 g_Norm: .real 0
