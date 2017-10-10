 .global RealMatrixMult
 .global ComplexMatrixMult
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

 .text

RealMatrixMult:
	move #0, a0.l
	move g_nRows1, a0.s
	move (a0), r0.l; R0 <- nRows1
	move g_nColumns1, a0.s
	move (a0), r2.l; R2 <- nColumns1 == nRows2
	move g_nColumns2, a0.s
	move (a0), r4.l; R4 <- nColumns2
	move g_Factor, a0.s
	move (a0), r18.l; R6 <- Factor
	move g_pReal1, a0.s
	move g_pReal2, a1.s
	move g_pReal3, a2.s
	; сохраняем исходные значения А0 и А1
	move a0, r14.s
	move a1, r16.s
	; задаём шаг адреса
	move r4.s, i1.s
	move #0, r20.l; счёткик строк
	RMMStartRowsLoop:
		move #0, r22.l; счёткик столбцов
		RMMStartColumnsLoop:
			; выставляем начальные адреса для цикла
			mpuu r20.s, r2.s, r24.l
			add r14.s, r24.s
			move r24.s, a0
			add r16.s, r22.s, r24.s
			move r24.s, a1
			move #0, r6.l; зануляем текущую сумму
			do r2, RMMEndInnerLoop
				move (a0)+, r8.l; копируем элемент первой матрицы в r8
				move (a1)+i1, r10.l; копируем элемент второй матрицы в r10
				fmpy r8.l, r10.l, r12.l; произведение
			RMMEndInnerLoop:
				fadd r12.l, r6.l; прибавляем его к сумме
			; кладём результат в матрицу, домножая на factor
			fmpy r18.l, r6.l
			move r6.l, (a2)+
			; инкрементим счётчик столбцов, проверяем условие выхода
			incl r22.l, r22.l
			cmpl r22.l, r4.l
			j.ne RMMStartColumnsLoop
		; инкрементим счётчик строк, проверяем условие выхода
		incl r20.l, r20.l
		cmpl r20.l, r0.l
		j.ne RMMStartRowsLoop
	stop


ComplexMatrixMult:
	move #0, a0.l
	move g_nRows1, a0.s
	move (a0), r0.l; nRows1
	move g_nColumns1, a0.s
	move (a0), r2.l; nColumns1 == nRows2
	move g_nColumns2, a0.s
	move (a0), r4.l; nColumns2
	move g_Factor, a0.s
	move (a0), r18.l; Factor
	move g_pReal1, a0.s
	move g_pReal2, a1.s
	move g_pReal3, a2.s
	move g_pImag1, a3.s
	move g_pImag2, a4.s
	move g_pImag3, a5.s
	; сохраняем исходные значения a0 и a1
	move a0, r14.s
	move a1, r16.s
	move a3, r26.s
	move a4, r28.s
	; задаём шаг адреса
	move r4.s, i1.s
	move r4.s, i4.s
	move #0, r20.l; счёткик строк
	CMMStartRowsLoop:
		move #0, r22.l; счёткик столбцов
		CMMStartColumnsLoop:
			; выставляем начальные адреса для цикла
			; re A
			mpuu r20.s, r2.s, r24.l; отступ в первой матрице = количество столбцов * счётчик строк
			add r14.s, r24.s, r30.s; прибавляем его к сохранённому в r14 значению a0
			move r30.s, a0; выставлям а0
			; im A
			add r26.s, r30.s; прибавляем отступ к сохранённому в r26 значению а3
			move r30.s, a3; выставляем а3
			; re B
			add r16.s, r22.s, r24.s; отступ во второй матрице = счётчик столбцов. прибавляем его к сохранённому в r16 значению а1
			move r24.s, a1
			; im B
			add r28.s, r22.s, r24.s; аналогично
			move r24.s, a4
			; зануляем текущую сумму
			move #0, r6.l
			move #0, r7.l
			do r2, CMMEndInnerLoop
				move (a0)+, r8.l; копируем элемент re первой матрицы в r8
				move (a1)+i1, r10.l; копируем элемент re второй матрицы в r10
				move (a3)+, r24.l; копируем элемент im первой матрицы в r24
				move (a4)+i4, r30.l; копируем элемент im второй матрицы в r30
				; re = re1*re2 - im1*im2
				fmpy r8.l, r10.l, r12.l
				fadd r12.l, r6.l
				fmpy r24.l, r30.l, r12.l
				fsub r12.l, r6.l
				; im = re1*im2 + re2*im1
				fmpy r8.l, r30.l, r30.l
				fadd r30.l, r7.l
				fmpy r10.l, r24.l, r24.l
			CMMEndInnerLoop:
				fadd r24.l, r7.l
			; кладём результат в матрицу, домножая на factor
			fmpy r18.l, r6.l
			move r6.l, (a2)+
			fmpy r18.l, r7.l, r6.l
			move r6.l, (a5)+
			; инкрементим счётчик столбцов, проверяем условие выхода
			incl r22.l, r22.l
			cmpl r22.l, r4.l
			j.ne CMMStartColumnsLoop
		; инкрементим счётчик строк, проверяем условие выхода
		incl r20.l, r20.l
		cmpl r20.l, r0.l
		j.ne CMMStartRowsLoop
	stop

 .data
 g_pReal2: .space 128, 0; 19208, 0
 g_pImag2: .space 128, 0; 19208, 0

 g_pReal1: .space 64, 0; 9604, 0
 g_pImag1: .space 64, 0; 9604, 0
 g_pReal3: .space 64, 0; 9604, 0
 g_pImag3: .space 128, 0; 9604, 0
 g_Factor: .real 0
 g_nRows1: .word 0
 g_nColumns1: .word 0
 g_nColumns2: .word 0



 g_FactorIm: .real 0
