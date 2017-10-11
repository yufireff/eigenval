.global ComplexMatrixInv
.global g_pRealDoubleSize
.global g_pImagDoubleSize

.text

ComplexMatrixInv:
	move #0, a0.l
	move g_nRows1, a0.s
	move (a0), r0.l ; r0.s = N
	mpuu #2, r0.s, r2.l
	move r2.s, r1.s ; r1.s = 2N
	; нельзя сравнивать счётчики цикла с r0.l, так как r0.l = N | (2N << 16)
; занятые регистры ri.l
; x 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	move #0, a0.l
	move g_pReal1, a0.s
	move #0, a1.l
	move g_pImag1, a1.s
	move #0, a2.l
	move g_pReal2, a2.s
	move #0, a3.l
	move g_pImag2, a3.s

	; запоминаем начальные адреса
	move a2, r2.s ; r2.s = a2
	move a3, r3.s ; r3.s = a3
; x 1 x 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	move #0, r4.l ; литерал 0.0f
	move 0x3f800000, r6.l ; литерал 1.0f
; x 1 x 3 x 5 x 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	; записываем присоединённую матрицу
	move #0, r8.s ; счётчик строк
	CMMStartAssociatedMatrixRowsLoop:
		move #0, r9.s ; счётчик столбцов
; x 1 x 3 x 5 x 7 x 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
		CMMStartAssociatedMatrixColumnsLoop1: ; часть матрицы с данными
			move (a0)+, r10.l
			move r10.l, (a2)+
			move (a1)+, r10.l
			move r10.l, (a3)+
			inc r9.s, r9.s
			cmp r9.s, r0.s
			j.ne CMMStartAssociatedMatrixColumnsLoop1

		move #0, r9.s
		CMMStartAssociatedMatrixColumnsLoop2: ; часть единичной матрицы
			cmp r8.s, r9.s
			move.eq r6.l, (a2)+
			move.ne r4.l, (a2)+
			move r4.l, (a3)+
			inc r9.s, r9.s
			cmp r9.s, r0.s
			j.ne CMMStartAssociatedMatrixColumnsLoop2

		inc r8.s, r8.s
		cmp r8.s, r0.s
		j.ne CMMStartAssociatedMatrixRowsLoop
; x 1 x 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	; возвращаем адреса на начальные значения
	move r2.s, a2
	move r3.s, a3
	move #1, i2.s
	move #1, i3.s


	; приводим присоединенную матрицу к верхнетреугольному виду
	move #0, r8.s ; счётчик строк i
	sub #1, r0.s, r7.s ; r7.s = N - 2
	move 0x40000000, r5.l ; литерал 2.0f
	CMMStartLowerLeftMainLoop:
		; запоминаем диагональный элемент
		mpuu r8.s, r1.s, r4.l ; r4.s = 2N * i
		add r4.s, r8.s, r6.s ; r6.s = 2N+i + i
; x 1 x 3 x x x 7 x 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 x 31
		add r2.s, r6.s, r10.s
		move r10.s, a2
		add r3.s, r6.s, r10.s
		move r10.s, a3
		move (a2), r10.l ; r10 = re_diag
		move (a3), r12.l ; r12 = im_diag

		; вычисляем норму диагонального элемента
		fmpy r10.l, r10.l, r1.l
		fmpy r12.l, r12.l, r3.l
		fadd r3.l, r1.l ; r1 = re^2 diag + im^2 diag
; x x x 3 x x x 7 x 9 x 11 x 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 x 31
		fin r1.l, r3.l ; r1.l - первое приближение 1 / norm(diag)
		; вычисляем второе приближение по формуле Ньютона-Рафсона y_{i+1} = y_i*(2 - x*y_i)
		; если точности не хватит, можно взять несколько итераций в цикле
		do 2, CMMEndNormDiag
			fmpy r1.l, r3.l, r7.l
			fsub r7.l, r5.l, r7.l
			CMMEndNormDiag:
			fmpy r7.l, r3.l
; x x x x x x x 7 x 9 x 11 x 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 x 31
		; запоминаем адреса начала строки (start + 2N*i)
		add r2.s, r4.s, r14.s
		add r3.s, r4.s, r15.s
; x x x x x x x 7 x 9 x 11 x 13 x 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 x 31

		; делаем диагональный элемент равным 1
		move r8.s, r16.s
; x x x x x x x 7 x 9 x 11 x 13 x 15 x 17 18 19 20 21 22 23 24 25 26 27 28 29 x 31
		CMMStartDiagonalOneLoop:
			; делим текущий элемент на диагональный
			move (a2), r18.l ; re_i
			move (a3), r20.l ; im_i
; x x x x x x x 7 x 9 x 11 x 13 x 15 x 17 x 19 x 21 22 23 24 25 26 27 28 29 x 31
			fmpy r18.l, r10.l, r7.l ; re i * re diag
			fmpy r20.l, r12.l, r9.l ; im_i * im_diag
			fadd r9.l, r7.l ; r7 = re_i * re_diag - im_i * im_diag
			fmpy r7.l, r3.l, r22.l; делим действительную часть на норму
			move r22.l, (a2)+ ; записываем результат

			fmpy r18.l, r12.l, r7.l ; re i * im diag
			fmpy r20.l, r10.l, r9.l ; im_i * re_diag
			fsub r7.l, r9.l ; r9 = re_diag * im_i - im_diag * re_i
			fmpy r9.l, r3.l, r22.l; делим мнимую часть на норму
			move r22.l, (a3)+ ; записываем результат

			inc r16.s, r16.s
			cmp r16.s, r1.s
			j.ne CMMStartDiagonalOneLoop

		; проверяем, не последняя ли это строка. если да, выходим досрочно
		cmp r8.s, r7.s
		j.eq CMMEndLowerLeftMainLoop

		; обнуляем остаток текущего столбца
		; для этого вычитаем из всех следующих строк текущую с коэффициентом
		; введём 2 адреса:
		; а0, а1 - для вычитаемой строки. начинается с base + 2N * i
		; a2, a3 - для текущей обнуляемой. начинается с base + 2N * j, j = [i+1, N]
		; начало вычитаемой строки сохранено в регистрах r14, r15
		; начинаем с первого ненулевого элемента
		inc r8.s, r22.s ; начинаем цикл с i + 1
; x x x x x x x 7 x 9 x 11 x 13 x 15 x 17 x 19 x 21 x 23 24 25 26 27 28 29 30 31
		CMMStartZeroingColumnOuterLoop1:
			add r14.s, r8.s, r24.s ; 2N*i + i
			move r24.s, a0
			add r15.s, r8.s, r24.s
			move r24.s, a1
			; задаём а2 и а3
			mpuu r22.s, r1.s, r24.l ; 2N*j
			add r8.s, r24.s ; 2N*j + i
			add r2.s, r24.s, r26.s
			move r26.s, a2.s
			add r3.s, r24.s, r26.s
			move r26.s, a3.s
			; запоминаем узловой элемент, на который умножается строка
			move (a2), r24.l ; двойной move для чтения в нечётные регистры - чётные заняты...
			move r24.l, r7.l ; re 1
			move (a3), r24.l
			move r24.l, r9.l ; im 1
; x x x x x x x x x x x 11 x 13 x 15 x 17 x 19 x 21 x 23 24 25 x 27 28 29 30 31
			; цикл по строке
			move r8.s, r24.s
; x x x x x x x x x x x 11 x 13 x 15 x 17 x 19 x 21 x 23 x 25 x 27 28 29 30 31
			CMMStartZeroingColumnInnerLoop1:
				; элемент вычитаемой строки
				move (a0)+, r26.l ; не хватает чётных регистров, поэтому одно из чисел читаем в 2 приёма
				move r26.l, r11.l ; re 2
				move (a1)+, r30.l ; im 2
				; текущий элемент строки
				move (a2), r26.l ; re 3
				move (a3), r28.l ; im 3
; x x x x x x x x x x x x x 13 x 15 x 17 x 19 x 21 x 23 x 25 x 27 x 29 30 x
				; произведение вычитаемого элемента на узловой
				fmpy r7.l, r11.l, r19.l ; re1 * re2
				fmpy r9.l, r30.l, r21.l ; im1 * im2
				fmpy r7.l, r30.l, r23.l ; re1 * im2
				fmpy r11.l, r9.l, r25.l ; re2 * im1
				fsub r21.l, r19.l ; re1 * re2 - im1 * im2
				fadd r25.l, r23.l ; re1 * im2 + re2 * im1
				; вычитаем его из текущего элемента
				fsub r19.l, r26.l
				fsub r23.l, r28.l
				; запоминаем результат
				move r26.l, (a2)+
				move r28.l, (a3)+

				inc r24.s, r24.s ; выход из цикла вычитания строки
				cmp r24.s, r1.s
				j.ne CMMStartZeroingColumnInnerLoop1
			inc r22.s, r22.s ; выход из цикла обнуления столбца
			cmp r22.s, r0.s
			j.ne CMMStartZeroingColumnOuterLoop1
		inc r8.s, r8.s
		cmp r8.s, r0.s
		j.ne CMMStartLowerLeftMainLoop
	CMMEndLowerLeftMainLoop:
; x 1 x 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	; обнуляем верхний треугольник
	; a0, a1 - адрес вычитаемой строки
	; a2, a3 - адрес текущей строки
	; занулять начинам с последнего столбца, i = N-1..1. строки от 0 до i-1
	; по строке идём от i до 2N
	dec r0.s, r4.s ; диагональный индекс i = N - 1
	CMMStartUpperRightMainLoop:
		dec r4.s, r5.s ; r5 = i - 1
		; задаём а0 и а1
		mpuu r4.s, r1.s, r6.l ; r6.s = i * 2N
		add r6.s, r4.s, r7.s ; r7 = i*2N + i
		add r2.s, r7.s, r8.s ; r8 = g_pReal2 + i*2N + i
		add r3.s, r7.s, r9.s ; r9 = g_pImag2 + i*2N + i
; x 1 x 3 x 5 x 7 x 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

		move #0, r10.s ; номер текущей строки j = 0..i-1
		CMMStartZeroingColumnOuterLoop2:
			move r8.s, a0.s
			move r9.s, a1.s
			; текущую строку начинаем с j*2N + i
			mpuu r10.s, r1.s, r12.l ; r12.s = j * 2N
			add r12.s, r4.s, r11.s; r11.s = j*2N + i
; x 1 x 3 x 5 x 7 x 9 x 11 x 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
			add r2.s, r11.s, r14.s
			move r14.s, a2 ; a2 = g_pReal2 + j*2N + i
			add r3.s, r11.s, r14.s
			move r14.s, a3 ; a3 = g_pImag2 + j*2N + i
			; запоминаем узловой элемент
			move (a2), r14.l ; re 1
			move (a3), r16.l ; im 1
; x 1 x 3 x 5 x 7 x 9 x 11 x 13 x 15 x 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
			move r4.s, r18.s ; счётчик по строке k = i..2N
			CMMStartZeroingColumnInnerLoop2:
				; запоминаем элемент вычитаемой строки
				move (a0)+, r20.l ; re 2
				move (a1)+, r22.l ; im 2
; x 1 x 3 x 5 x 7 x 9 x 11 x 13 x 15 x 17 x 19 x 21 x 23 24 25 26 27 28 29 30 31
				; умножаем его на узловой элемент
				fmpy r14.l, r20.l, r1.l ; re1 * re2
				fmpy r16.l, r22.l, r3.l ; im1 * im2
				fmpy r14.l, r22.l, r5.l ; re1 * im2
				fmpy r20.l, r16.l, r7.l ; re2 * im1
				fsub r3.l, r1.l ; r1.l = re1*re2 - im1*im2
				fadd r7.l, r5.l ; r5.l = re1 * im2 + re2 * im1
				; вычитаем его из текущего элемента
				move (a2), r24.l ; re 3
				move (a3), r26.l ; im 3
				fsub r1.l, r24.l ; re3 - ...
				fsub r5.l, r26.l ; im3 - ...
				; запоминаем результат
				move r24.l, (a2)+
				move r26.l, (a3)+
				inc r18.s, r18.s ; выход из цикла вычитания строки
				cmp r18.s, r1.s
				j.ne CMMStartZeroingColumnInnerLoop2
			inc r10.s, r10.s ; выход из цикла зануления столбца
			cmp r10.s, r4.s
			j.ne CMMStartZeroingColumnOuterLoop2
		dec r4.s, r4.s ; выход из цикла по диагональным элементам
		cmp #0, r4.s
		j.ne CMMStartUpperRightMainLoop
; x 1 x 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	; обратная матрица найдена. копируем результат в буфер 1
	move g_pReal1, a0.s
	move g_pImag1, a1.s
	add r2.s, r0.s, r2.s
	add r3.s, r0.s, r3.s
	move #0, r4.s ; i
	CMMStartCopyRowLoop:
		move r2.s, a2.s
		move r3.s, a3.s
		move #0, r6.s
		CMMStartCopyColumnLoop:
			move (a2)+, r8.l
			move r8.l, (a0)+
			move (a3)+, r8.l
			move r8.l, (a1)+
			inc r6.s, r6.s
			cmp r6.s, r0.s
			j.ne CMMStartCopyColumnLoop
		add r2.s, r1.s, r2.s
		add r3.s, r1.s, r3.s
		inc r4.s, r4.s
		cmp r4.s, r0.s
		j.ne CMMStartCopyRowLoop
	stop
