 .global RealMatrixMultAtBA
 .global MatrixAminusUSUt
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

; умножение действительных матриц A' * B * A
RealMatrixMultAtBA:
	move g_Buffer1, a0.s
	move g_Buffer3d, a1.s
	move g_Buffer5, a2.s
	move g_Buffer2, a3.s ; используем буфер Imag1 для промежуточных вычислений
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

; умножение комплексных матриц A - U * S * U', где S - диагональная действительная
MatrixAminusUSUt:
	; U расположена в g_Buffer1 и g_Buffer2
	move g_Buffer1, a0.s
	move g_Buffer2, a1.s
	; S расположена в m_pDiag
	move g_pDiag, a2.s
	; A расположена в g_Buffer3d и g_Buffer4d
	; промежуточное произведение считаем в g_Buffer5 и g_Buffer6
	move g_Buffer5, a5.s
	move g_Buffer6, a6.s
	; все матрицы квадратные, размера g_Size
	move g_Size, a7.s
	move (a7), r8.l ; r8.s = N
	; квадрат нормы считаем в g_Norm

	; начальное состояние адресов
	move a0.s, r0.s
	move a1.s, r1.s
	move a2.s, r2.s
	move a3.s, r3.s
	move a4.s, r4.s
	move a5.s, r5.s
	move a6.s, r6.s
; занятые регистры ri.l
; x 1 x 3 x 5 x 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	; сначала считаем произведение U * S, где U - квадратная комплексная, S - диагональная действительная
	; US[i, j] = U[i, j] * S[j]
	move #0, r7.s ; i
	MAmUSUtStartDiagMultRowsLoop:
		move #0, r9.s ; j
		mpuu r7.s, r8.s, r10.l
		add r10.s, r0.s, r12.s
		move r12.s, a0.s
		add r10.s, r1.s, r12.s
		move r12.s, a1.s
		move r2.s, a2.s
		MAmUSUtStartDiagMultColumnsLoop:
			move (a0)+, r10.l ; re U[i, j]
			move (a1)+, r12.l ; im U[i, j]
			move (a2)+, r14.l ; S[j]
			; перемножаем
			fmpy r10.l, r14.l, r16.l
			fmpy r12.l, r14.l, r18.l
			; записываем результат
			move r16.l, (a5)+
			move r18.l, (a6)+
			inc r9.s, r9.s ; выход из цикла по столбцам
			cmp r9.s, r8.s
			j.ne MAmUSUtStartDiagMultColumnsLoop
		inc r7.s, r7.s ; выход из цикла по строкам
		cmp r7.s, r8.s
		j.ne MAmUSUtStartDiagMultRowsLoop

	move r0.s, a0.s
	move r1.s, a1.s
	move r5.s, a5.s
	move r6.s, a6.s
	move g_Buffer7, a3.s
	move g_Buffer8, a4.s
; x 1 x 3 x 5 x 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	; результат умножаем на U'
	; в обеих матрицах идём по строкам
	; результат записываем в g_Buffer7, g_Buffer8
	move #0, r7.s ; i
	MAmUSUtStartTranspMultRowsLoop:
		move #0, r9.s ; j
		MAmUSUtStartTranspMultColumnsLoop:
			; начальное значение a5, a6 = base + i * N
			mpuu r8.s, r7.s, r10.l
			add r5.s, r10.s, r12.s
			move r12.s, a5.s
			add r6.s, r10.s, r12.s
			move r12.s, a6.s
			; начальное значение a0, a1 = base + j * N
			mpuu r8.s, r9.s, r10.l
			add r0.s, r10.s, r12.s
			move r12.s, a0.s
			add r1.s, r10.s, r12.s
			move r12.s, a1.s
			; частичная сумма
			move #0, r10.l ; re
			move #0, r12.l ; im
			do r8.s, MAmUSUtEndInnerLoop
				move (a5)+, r14.l ; re 1
				move (a6)+, r16.l ; im 1
				move (a0)+, r18.l ; re 2
				move (a1)+, r20.l ; im 2
				fmpy r14.l, r18.l, r1.l ; re1 * re2
				fmpy r16.l, r20.l, r3.l ; im1 * im2
				fmpy r14.l, r20.l, r5.l ; re1 * im2
				fmpy r18.l, r16.l, r7.l ; re2 * im1
				fadd r3.l, r1.l, r9.l ; re1*re2 + im1*im2 (знак im2 изменён из-за транспонирования)
				fsub r5.l, r7.l, r11.l ; re2*im1 - re1*im2  (знак im2 изменён из-за транспонирования)
				fadd r9.l, r10.l ; re += ...
			MAmUSUtEndInnerLoop:
				fadd r11.l, r12.l; im += ...
			; запоминаем значения
			move r10.l, (a3)+
			move r12.l, (a4)+
			inc r9.s, r9.s ; выход из цикла по столбцам
			cmp r9.s, r8.s
			j.ne MAmUSUtStartTranspMultColumnsLoop
		inc r7.s, r7.s ; выход из цикла по строкам
		cmp r7.s, r8.s
		j.ne MAmUSUtStartTranspMultRowsLoop

	; вычитаем полученную матрицу из А, результат кладём в g_Buffer5, g_Buffer6
	; попутно считаем квадрат нормы матрицы
	; A
	move g_Buffer3d, a0.s
	move g_Buffer4d, a1.s
	; U * S * U'
	move g_Buffer7, a3.s
	move g_Buffer8, a4.s
	; A - U * S * U'
	move g_Buffer5, a5.s
	move g_Buffer6, a6.s

	move #0, r30.l ; норма
	mpuu r8.s, r8.s, r10.l ; N*N
	do r10.s, MAmUSUtEndSubtractLoop
		move (a0)+, r12.l ; re 1
		move (a1)+, r14.l ; im 1
		move (a3)+, r16.l ; re 2
		move (a4)+, r18.l ; im 2
		fsub r16.l, r12.l, r20.l ; re 1 - re 2
		fsub r18.l, r14.l, r22.l ; im 1 - im 2
		; записываем результат
		move r20.l, (a5)+
		move r22.l, (a6)+
		; считаем норму
		fmpy r20.l, r20.l, r1.l
		fmpy r22.l, r22.l, r3.l
		fadd r1.l, r30.l
	MAmUSUtEndSubtractLoop:
		fadd r3.l, r30.l

	; записываем норму в g_Norm
	move g_Norm, a7.s
	move r30.l, (a7)
	stop
