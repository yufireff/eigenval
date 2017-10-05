.global ComplexMatrixInv
.global g_pRealDoubleSize
.global g_pImagDoubleSize

.text

ComplexMatrixInv:
	move #0, a0.l
	move g_nRows1, a0.s
	move (a0), r0.l ; N
	mpuu #2, r0.s, r16.l ; 2N

	move #0, a0.l
	move g_pReal1, a0.s
	move #0, a1.l
	move g_pImag1, a1.s
	move #0, a2.l
	move g_pReal2, a2.s
	move #0, a3.l
	move g_pImag2, a3.s

    mpuu r0.s, r16.s, r2.l ; N * 2N

	; запоминаем начальные адреса
	move a0, r3.s
	move a1, r4.s
	move a2, r5.s
	move a3, r6.s

	move #0, r10.l; eye
	move 0x3f800000, r12.l; 1.0f
	move #0, r7.l; счётчик строк
	; записываем присоединённую матрицу
	CMMStartAssociatedMatrixRowsLoop:
		move #0, r8.l; счётчик столбцов
		CMMStartAssociatedMatrixColumnsLoop1: ; часть единичной матрицы
			cmpl r7.l, r8.l
			move.eq r12.l, (a2)+
			move.ne r10.l, (a2)+
			move r10.l, (a3)+
			incl r8.l, r8.l
			cmpl r8.l, r0.l
			j.ne CMMStartAssociatedMatrixColumnsLoop1
		CMMStartAssociatedMatrixColumnsLoop2: ; часть матрицы с данными
			move (a0)+, r14.l
			move r14.l, (a2)+
			move (a1)+, r14.l
			move r14.l, (a3)+
			incl r8.l, r8.l
			cmpl r8.l, r16.l
			j.ne CMMStartAssociatedMatrixColumnsLoop2
		incl r7.l, r7.l
		cmpl r7.l, r0.l
		j.ne CMMStartAssociatedMatrixRowsLoop
	stop