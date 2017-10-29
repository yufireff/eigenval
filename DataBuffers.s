 .global g_Buffer1
 .global g_Buffer2
 .global g_Buffer3d
 .global g_Buffer4d
 .global g_Buffer5
 .global g_Buffer6
 .global g_Buffer7
 .global g_Buffer8
 .global g_pDiag
 .global g_Factor
 .global g_FactorIm
 .global g_nRows1
 .global g_nColumns1
 .global g_nRows2
 .global g_nColumns2
 .global g_Size
 .global g_Norm
 .global g_p
 .global g_q
 .global g_ptr0

 .data
 g_Buffer1: .space 49*49*4, 0 ; T
 g_Buffer2: .space 49*49*4, 0 ; U (Q)
 g_Buffer3d: .space 2*49*49*4, 0
 g_Buffer4d: .space 2*49*49*4, 0
 g_Buffer5: .space 49*49*4, 0 ; D (S)
 g_Buffer6: .space 49*49*4, 0 ; Qm
 g_Buffer7: .space 49*49*4, 0 ; Rm
 g_Buffer8: .space 49*49*4, 0
 g_pDiag:  .space 49*4, 0
 ; для отладки
 ;g_Buffer1: .space 4*4*4, 0
 ;g_Buffer2: .space 4*4*4, 0
 ;g_Buffer3d: .space 4*4*4, 0
 ;g_Buffer4d: .space 4*4*4, 0
 ;g_Buffer5: .space 4*4*4, 0
 ;g_Buffer6: .space 4*4*4, 0
 ;g_pDiag:  .space 4*4, 0
 ;g_Buffer7: .space 4*4*4, 0
 ;g_Buffer8: .space 4*4*4, 0

 g_Factor: .real 0
 g_nRows1: .word 0
 g_nColumns1: .word 0
 g_nRows2: .word 0
 g_nColumns2: .word 0
 g_FactorIm: .real 0
 g_Size: .word 0
 g_p: .word 0
 g_q: .word 0
 g_Norm: .real 0
 g_ptr0: .word 0
