.global ComplexMatrixInv
.global g_pRealDoubleSize
.global g_pImagDoubleSize

.text

ComplexMatrixInv:
	move #0, a0.l
	move g_nRows1, a0.s
	move (a0), r0.l ; N
	mpuu #2, r0.s, r4.l ; 2N

	move #0, a0.l
	move g_pReal1, a0.s
	move #0, a1.l
	move g_pImag1, a1.s
	move #0, a2.l
	move g_pReal2, a2.s
	move #0, a3.l
	move g_pImag2, a3.s

    mpuu r0.s, r4.s, r2.l ; N * 2N

	; запоминаем начальные адреса
	move a2, r5.l
	move a3, r6.l

	move #0, r8.l; eye
	move 0x3f800000, r10.l; 1.0f

	; зан€ты регистры r0, r2, r4, r5, r6, r8, r10

	move #0, r1.l; счЄтчик строк
	; записываем присоединЄнную матрицу
	CMMStartAssociatedMatrixRowsLoop:
		move #0, r3.l; счЄтчик столбцов
		CMMStartAssociatedMatrixColumnsLoop1: ; часть матрицы с данными
			move (a0)+, r14.l
			move r14.l, (a2)+
			move (a1)+, r14.l
			move r14.l, (a3)+
			incl r3.l, r3.l
			cmpl r3.l, r0.l
			j.ne CMMStartAssociatedMatrixColumnsLoop1
			move #0, r3.l
		CMMStartAssociatedMatrixColumnsLoop2: ; часть единичной матрицы
			cmpl r1.l, r3.l
			move.eq r10.l, (a2)+
			move.ne r8.l, (a2)+
			move r8.l, (a3)+
			incl r3.l, r3.l
			cmpl r3.l, r0.l
			j.ne CMMStartAssociatedMatrixColumnsLoop2

		incl r1.l, r1.l
		cmpl r1.l, r0.l
		j.ne CMMStartAssociatedMatrixRowsLoop

	; возвращаем адреса на начальные значени€
	move r5.l, a2
	move r6.l, a3

	; приводим присоединенную матрицу к верхнетреугольному виду
	move #0, r1.l; счЄтчик строк
	subl #2, r0.l, r25.l ; r25 = N - 2
	CMMStartLowerLeftMainLoop:
		; запоминаем диагональный элемент
		mpuu r1.s, r4.s, r3.l ; r3 = 2N * i
		nop ; умножение происходит в 2 такта
		addl r1.l, r3.l, r7.l ; r7 = 2N*i + i
		move r3.s, i2.s
		move r3.s, i3.s
		move (a2 + i2), r12.l ; re diag
		move (a3 + i3), r14.l ; im diag
		; делаем диагональный элемент равным 1
		; вычисл€ем норму диагонального элемента
		fmpy r12.l, r12.l, r13.l ; re^2 diag
		fmpy r14.l, r14.l, r15.l ; im^2 diag. складывать будем позже, чтобы не вставл€ть nop
		move #0, r9.l ; обнул€ем индекс цикла
		; запоминаем адреса начала строки
		addl r5.l, r3.l, r11.l
		move r11.l, a2
		addl r6.l, r3.l, r27.l
		fadd r15.l, r13.l ; r13 = re^2 diag + im^2 diag. сложение перенесено сюда в цел€х оптимизации
		move r27.l, a3
		move #1, i2.s
		move #1, i3.s
		fin r13.l, r13.l ; r13 = 1 / (re^2 diag + im^2 diag). деление перенесено также в цел€х оптимизации
		CMMStartDiagonalOneLoop:
			; делим текущий элемент на диагональный
			move (a2), r16.l ; re
			move (a3), r18.l ; im
			fmpy r12.l, r16.l, r15.l ; re_i * re_diag
			fmpy r14.l, r18.l, r17.l ; im_i * im_diag
			;nop
			fadd r17.l, r15.l ; r15 = re_i * re_diag + im_i * im_diag
			; чтобы было меньше nop, считает другие произведени€
			fmpy r12.l, r18.l, r19.l ; re diag * im i
			fmpy r16.l, r14.l, r21.l ; re i * im diag
			fmpy r15.l, r13.l, r20.l ; делим действительную часть на норму
			fsub r19.l, r21.l, r23.l ; числитель мнимой части
			move r20.l, (a2)+ ; пока вычитаетс€ числитель мнимой части, записываем действительную
			fmpy r23.l, r13.l, r20.l ; делим мнимую часть на норму
			incl r9.l, r9.l ; пока делитс€ мнима€ часть, инкрементим счЄтчик цикла
			move r22.l, (a3)+ ; записываем мнимую часть в матрицу

			; провер€ем, не последн€€ ли это строка. если да, выходим досрочно
			cmpl r1.l, r25.l
			j.eq CMMEndLowerLeftMainLoop

			; обнул€ем остаток текущего столбца
			; дл€ этого вычитаем из всех следующих строк текущую с коэффициентом
			; введЄм 2 адреса:
			; а0, а1 - дл€ вычитаемой строки. начинаетс€ с base + 2N * i
			; a2, a3 - дл€ текущей обнул€емой. начинаетс€ с base + 2N * j, j = [i+1, N]
			; начало вычитаемой строки сохранено в регистрах r11, r27
			move r11.l, a0
			move r27, a1
			
			
			cmpl r9.l, r4.l ; провер€ем условие выхода из цикла
			j.ne CMMStartDiagonalOneLoop
		incl r1.l, r1.l
		cmpl r1.l, r0.l
		j.ne CMMStartLowerLeftMainLoop
	CMMEndLowerLeftMainLoop:
	
	stop