 .global RealMatrixMultLeftTransp
 .global ComplexMatrixMultLeftTransp
 .global ComplexMatrixMultRightTransp
 .global g_Buffer1
 .global g_Buffer2
 .global g_Buffer3d
 .global g_Buffer4d
 .global g_Buffer5
 .global g_Buffer6
 .global g_Factor
 .global g_FactorIm
 .global g_nRows1
 .global g_nColumns1
 .global g_nRows2
 .global g_nColumns2

 .text

; умножение действительных матриц factor * A' * B
RealMatrixMultLeftTransp:
	move g_Buffer1, a0.s
	move g_Buffer3d, a1.s
	move g_Buffer5, a2.s
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
			move #0, r18.l ; частичная сумма
			do r4.s, RMMLTEndInnerLoop
				move (a0)+i0, r12.l
				move (a1)+i1, r14.l
				fmpy r12.l, r14.l, r16.l
			RMMLTEndInnerLoop:
				fadd r16.l, r18.l
				fmpy r10.l, r18.l, r16.l ; умножаем на factor
				move r16.l, (a2)+
			inc r9.s, r9.s ; условие выхода из цикла по столбцам
			cmp r9.s, r8.s
			j.ne RMMLTStartColumnsLoop
		inc r7.s, r7.s ; условие выхода из цикла по строкам
		cmp r7.s, r6.s
		j.ne RMMLTStartRowsLoop
	stop

; умножение комплексных матриц A' * B
ComplexMatrixMultLeftTransp:
	move g_Buffer1, a0.s
	move g_Buffer2, a1.s
	move g_Buffer3d, a2.s
	move g_Buffer4d, a3.s
	move g_Buffer5, a4.s
	move g_Buffer6, a5.s
	; начальные значения адресов
	move a0.s, r0.s
	move a1.s, r1.s
	move a2.s, r2.s
	move a3.s, r3.s
	move a4.s, r4.s
	move a5.s, r5.s

	move g_nRows1, a6.s
	move (a6), r6.l; r6.s = nRows1 == nRows2
	move g_nColumns1, a6.s
	move (a6), r8.l; r8.s = nColumns1 == nRows3
	move g_nColumns2, a6.s
	move (a6), r10.l; r10.s = nColumns2 == nColumns3
; занятые регистры ri.l
; x 1 x 3 x 5 x 7 x 9 x 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	; задаём шаг адреса
	; для а0 - а3 - количество столбцов в соответствующих матрицах
	; для а4, а5 - 1
	move r8.s, i0.s
	move r8.s, i1.s
	move r10.s, i2.s
	move r10.s, i3.s
	move #1, i4.s
	move #1, i5.s

	move #0, r9.s ; i, счётчик строк в результирующей матрице
	СMMLTStartRowsLoop: ; в обеих матрицах идём вдоль столбцов
		move #0, r11.s ; j, счётчик столбцов в результирующей матрице
		СMMLTStartColumnsLoop:
			; начальное положение адреса a0, а1 = base + i
			add r0.s, r9.s, r12.s
			move r12.s, a0.s
			add r1.s, r9.s, r12.s
			move r12.s, a1.s
			; начальное положение адреса a2, а3 = base + j
			add r2.s, r11.s, r12.s
			move r12.s, a2.s
			add r3.s, r11.s, r12.s
			move r12.s, a3.s
			move #0, r20.l ; частичная сумма, re
			move #0, r22.l ; частичная сумма, im
; x 1 x 3 x 5 x 7 x 9 x 11 12 13 14 15 16 17 x 19 20 21 22 23 24 25 26 27 28 29 30 31
			do r6.s, CMMLTEndInnerLoop
				move (a0)+i0, r12.l ; re 1
				move (a1)+i1, r14.l ; im 1
				move (a2)+i2, r16.l ; re 2
				move (a3)+i3, r18.l ; im 2
				fmpy r12.l, r16.l, r1.l ; re1 * re2
				fmpy r14.l, r18.l, r3.l ; im1 * im2
				fmpy r12.l, r18.l, r5.l ; re1 * im2
				fmpy r16.l, r14.l, r7.l ; re2 * im1
				fadd r3.l, r1.l, r9.l ; re1*re2 + im1*im2 (знак im1 изменён из-за транспонирования)
				fsub r7.l, r5.l, r11.l ; re1*im2 - re2*im1 (знак im1 изменён из-за транспонирования)
				fadd r9.l, r20.l ; re += ...
			CMMLTEndInnerLoop:
				fadd r11.l, r22.l; im += ...
			; запоминаем значения
			move r20.l, (a4)+
			move r22.l, (a5)+
			inc r11.s, r11.s ; условие выхода из цикла по столбцам
			cmp r11.s, r10.s
			j.ne СMMLTStartColumnsLoop
		inc r9.s, r9.s ; условие выхода из цикла по строкам
		cmp r9.s, r8.s
		j.ne СMMLTStartRowsLoop
	stop


; умножение комплексных матриц A * B'
ComplexMatrixMultRightTransp:
	move g_Buffer1, a0.s
	move g_Buffer2, a1.s
	move g_Buffer3d, a2.s
	move g_Buffer4d, a3.s
	move g_Buffer5, a4.s
	move g_Buffer6, a5.s
	; начальные значения адресов
	move a0.s, r0.s
	move a1.s, r1.s
	move a2.s, r2.s
	move a3.s, r3.s
	move a4.s, r4.s
	move a5.s, r5.s

	move g_nRows1, a6.s
	move (a6), r6.l; r6.s = nRows1 == nRows3
	move g_nColumns1, a6.s
	move (a6), r8.l; r8.s = nColumns1 == nColumns2
	move g_nRows2, a6.s
	move (a6), r10.l; r10.s = nRows2 == nColumns3
; занятые регистры ri.l
; x 1 x 3 x 5 x 7 x 9 x 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	move #0, r7.s ; i, счётчик строк в результирующей матрице
	СMMRTStartRowsLoop: ; в обеих матрицах идём вдоль строк
		move #0, r11.s ; j, счётчик столбцов в результирующей матрице
		СMMRTStartColumnsLoop:
			; начальное положение адреса a0, а1 = base + i*nColumns1
			mpuu r8.s, r7.s, r14.l
			add r0.s, r14.s, r12.s
			move r12.s, a0.s
			add r1.s, r14.s, r12.s
			move r12.s, a1.s
			; начальное положение адреса a2, а3 = base + j*nColumns2
			mpuu r8.s, r11.s, r14.l
			add r2.s, r14.s, r12.s
			move r12.s, a2.s
			add r3.s, r14.s, r12.s
			move r12.s, a3.s
			move #0, r20.l ; частичная сумма, re
			move #0, r22.l ; частичная сумма, im
; x 1 x 3 x 5 x 7 x 9 x 11 12 13 14 15 16 17 x 19 20 21 22 23 24 25 26 27 28 29 30 31
			do r8.s, CMMRTEndInnerLoop
				move (a0)+, r12.l ; re 1
				move (a1)+, r14.l ; im 1
				move (a2)+, r16.l ; re 2
				move (a3)+, r18.l ; im 2
				fmpy r12.l, r16.l, r1.l ; re1 * re2
				fmpy r14.l, r18.l, r3.l ; im1 * im2
				fmpy r12.l, r18.l, r5.l ; re1 * im2
				fmpy r16.l, r14.l, r7.l ; re2 * im1
				fadd r3.l, r1.l, r9.l ; re1*re2 + im1*im2 (знак im2 изменён из-за транспонирования)
				fsub r5.l, r7.l, r11.l ; re2*im1 - re1*im2  (знак im2 изменён из-за транспонирования)
				fadd r9.l, r20.l ; re += ...
			CMMRTEndInnerLoop:
				fadd r11.l, r22.l; im += ...
			; запоминаем значения
			move r20.l, (a4)+
			move r22.l, (a5)+
			inc r11.s, r11.s ; условие выхода из цикла по столбцам
			cmp r11.s, r10.s
			j.ne СMMRTStartColumnsLoop
		inc r7.s, r7.s ; условие выхода из цикла по строкам
		cmp r7.s, r6.s
		j.ne СMMRTStartRowsLoop
	stop
