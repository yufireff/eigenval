 .global RealMatrixMultAtBA
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

; умножение действительных матриц A' * B * A
RealMatrixMultAtBA:
	move g_pReal1, a0.s
	move g_pReal2, a1.s
	move g_pReal3, a2.s
	move g_pImag1, a3.s ; используем буфер Imag1 для промежуточных вычислений
	; начальные значеия адресов
	move a0.s, r0.s
	move a1.s, r1.s
	move a2.s, r2.s
	move a3.s, r3.s

	move g_Size, a4.l
    move (a4), r4.l
; занятые регистры ri.l
; x 1 x 3 x 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	; сначала считаем произведение A' * B
	; задаём шаг адреса
	; для а0 и а1 - размер матрицы
	; для а3 - 1
	move r4.s, i0.s
	move r4.s, i1.s
	move #1, i3.s

	move #0, r7.s ; i, счётчик строк в результирующей матрице
	RMMAtBAStartRowsLoop1: ; в обеих матрицах идём вдоль столбцов
		move #0, r9.s ; j, счётчик столбцов в результирующей матрице
		RMMAtBAStartColumnsLoop1:
			; начальное положение адреса a0 = base + i
			add r0.s, r7.s, r12.s
			move r12.s, a0.s
			; начальное положение адреса a1 = base + j
			add r1.s, r9.s, r12.s
			move r12.s, a1.s
			move #0, r18.l ; частичная сумма
			do r4.s, RMMAtBAEndInnerLoop1
				move (a0)+i0, r12.l
				move (a1)+i1, r14.l
				fmpy r12.l, r14.l, r16.l
			RMMAtBAEndInnerLoop1:
				fadd r16.l, r18.l
				move r18.l, (a3)+
			inc r9.s, r9.s ; условие выхода из цикла по столбцам
			cmp r9.s, r4.s
			j.ne RMMAtBAStartColumnsLoop1
		inc r7.s, r7.s ; условие выхода из цикла по строкам
		cmp r7.s, r4.s
		j.ne RMMAtBAStartRowsLoop1

	; теперь результат умножаем на А
	move #1, i2.s
	move #0, r5.s ; i, счётчик строк в результирующей матрице
	RMMAtBAStartRowsLoop2:
		move #0, r6.s ; j, счётчик столбцов в результирующей матрице
		RMMAtBAStartColumnsLoop2:
			; начальное положение адреса a3 = base + i* N
			mpuu r5.s, r4.s, r8.l
			add r3.s, r8.s
			move r8.s, a3.s
			; начальное положение адреса a0 = base + j
			add r0.s, r6.s, r8.s
			move r8.s, a0.s
			; частичная сумма
			move #0, r10.l
; x 1 x 3 x 5 x 7 x 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
			do r4, RMMAtBAEndInnerLoop2
				move (a3)+, r12.l
				move (a0)+i0, r14.l
				fmpy r12.l, r14.l, r3.l
			RMMAtBAEndInnerLoop2:
				fadd r3.l, r10.l
			; сохраняем результат
			move r10.l, (a2)+

			inc r6.s, r6.s ; условие выхода из цикла по столбцам
			cmp r6.s, r4.s
			j.ne RMMAtBAStartColumnsLoop2
		inc r5.s, r5.s ; условие выхода из цикла по строкам
		cmp r5.s, r4.s
		j.ne RMMAtBAStartRowsLoop2
	stop

