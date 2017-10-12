 .global RealMatrixMultLeftTransp
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
 .global g_nRows2
 .global g_nColumns2

 .text

; умножение матриц A' * B
RealMatrixMultLeftTransp:
	move g_pReal1, a0.s
	move g_pReal2, a1.s
	move g_pReal3, a2.s
	; начальные значеия адресов
	move a0.s, r0.s
	move a1.s, r1.s
	move a2.s, r2.s

	move g_nRows1, a3.s
	move (a3), r4.l; r4.s = nRows1 == nRows2
	move g_nColumns1, a3.s
	move (a3), r6.l; r6.s = nColumns1 == nRows3
	move g_nColumns2, a3.s
	move (a3), r8.l; r8.s = nColumns2 == nColumns3
	move g_Factor, a3.s
	move (a3), r10.l; r10.l <- Factor
; занятые регистры ri.l
; x 1 x 3 x 5 x 7 x 9 x 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	; задаём шаг адреса
	; для а0 и а1 - количество столбцов в соответствующих матрицах
	; для а2 - 1
	move r6.s, i0.s
	move r8.s, i1.s
	move #1, i2.s

	move #0, r7.s ; i, счётчик строк в результирующей матрице
	RMMLTStartRowsLoop: ; в обеих матрицах идём вдоль столбцов
		move #0, r9.s ; j, счётчик столбцов в результирующей матрице
		RMMLTStartColumnsLoop:
			; начальное положение адреса a0 = base + i
			add r0.s, r7.s, r12.s
			move r12.s, a0.s
			; начальное положение адреса a1 = base + j
			add r1.s, r9.s, r12.s
			move r12.s, a1.s
			move #0, r18.l
			do r4.s, CMMLTEndInnerLoop
				move (a0)+i0, r12.l
				move (a1)+i1, r14.l
				fmpy r12.l, r14.l, r16.l
			CMMLTEndInnerLoop:
				fadd r16.l, r18.l
				fmpy r10.l, r18.l, r16.l
				move r16.l, (a2)+
			inc r9.s, r9.s ; условие выхода из цикла по столбцам
			cmp r9.s, r8.s
			j.ne RMMLTStartColumnsLoop
		inc r7.s, r7.s ; условие выхода из цикла по строкам
		cmp r7.s, r6.s
		j.ne RMMLTStartRowsLoop
	stop
