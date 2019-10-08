 .global QRSymmReal

; вызов подпрограмм не производится. тело функций qr_step и giv включено в текст основной программы

 .text

QRSymmReal:
; T находится в g_Buffer1
; n - в g_Size
; tol - в g_Factor
; на выходе получаем матрицы U (Q, g_Buffer7), D (S, g_Buffer5)
	move g_Size, a0.s
	move (a0), r0.l ; N
	move g_Factor, a0.s
	move (a0), r2.l ; tol
	move g_Buffer1, a0.s ; T
	move g_Buffer7, a1.s ; U (Q)
	move g_Buffer5, a2.s ; D (S)

	move #0, r30.l ; литерал 0.0F
	move 0x3F000000, r31.l ; 0.5F
	move 0x3FC00000, r29.l ; 1.5F
	move 0x3F800000, r28.l ; 1.0F
	move 0x40000000, r27.l ; 2.0F
	move 0x34000000, r25.l ; FLOAT_EPS

	; U = eye(N,N,'single');
	move #0, r14.s
	QRSREyeOuterLoop0Begin:
		move r28.l, (a1)+
		do r0.s, QRSREyeInnerLoop0End
		QRSREyeInnerLoop0End:
			move r30.l, (a1)+
		inc r14.s, r14.s
		cmp r14.s, r0.s
		j.ne QRSREyeOuterLoop0Begin
	move r28.l, (a1)

	; p = 1;
	move #1, r4.s ; p
	; q = N;
	move r0.s, r5.s ; q
; занятые регистры ri.l
; x 1 x 3 x 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 x 26 x x x x x

	; while q ~= 1
	QRSRStartMainLoop:
		; fl = 0;
		move #0, r6.s ; f1
		; for i = q-1:-1:1
		dec r5.s, r7.s ; i
; x 1 x 3 x x x 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 x 26 x x x x x
		QRSRStartInnerLoop:
; здесь индексация MatLab, начинается с 1
;			if (abs(T(i+1,i)) <= tol*(abs(T(i,i)) + abs(T(i+1,i+1))))
; чтобы правильно перевести указанное выражение, запишем его на языке С
; if (f_abs(Tcopy.pData[i*n + i-1]) <=
;			tol * (f_abs(Tcopy.pData[(i-1)*n + i-1]) + f_abs(Tcopy.pData[i*n + i])))
			move g_Buffer1, a0.s ; T
			dec r7.s, r8.s ; i-1
; x 1 x 3 x x x 7 x 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 x 26 x 28 x x 31
			mpuu r7.s, r0.s, r10.l ; i*n
			mpuu r8.s, r0.s, r12.l ; (i-1)*n
; x 1 x 3 x x x 7 x 9 x 11 x 13 14 15 16 17 18 19 20 21 22 23 24 x 26 x x x x x
			add r10.s, r7.s, r11.s ; i*n + i
			move r11.s, i0.s
			move (a0+i0), r16.l ; T(i+1,i+1)	; Tcopy.pData[i*n + i]
			andl 0x7fffffff, r16.l, r1.l ; abs(T(i+1,i+1))

			dec r11.s, r13.s
			move r13.s, i0.s
			move (a0+i0), r16.l ; T(i+1,i)	; Tcopy.pData[i*n + i-1]
			andl 0x7fffffff, r16.l, r5.l ; abs(T(i+1,i))

			add r12.s, r8.s, r15.s
			move r15.s, i0.s
			move (a0+i0), r16.l ; T(i,i)	; Tcopy.pData[(i-1)*n + i-1]
			andl 0x7fffffff, r16.l, r3.l ; abs(T(i,i))

			fadd r3.l, r1.l, r7.l ; abs(T(i,i)) + abs(T(i+1,i+1))
			fmpy r2.l, r7.l, r9.l ; tol*(abs(T(i,i)) + abs(T(i+1,i+1)))
			fsub r9.l, r5.l, r11.l ; abs(T(i+1,i)) - tol*(abs(T(i,i)) + abs(T(i+1,i+1)))
			; проверяем знак получившегося числа. если старший бит = 1, число отрицательное
			andl 0x80000000, r11.l, r19.l; если результат > 0, то разность отрицательна, и изначальное условие false
			j.eq QRSRCond1False
				; T(i,i+1) = 0;	Tcopy.pData[(i-1)*n + i]
				inc r15.s, r17.s
				move r17.s, i0.s
				move r30.l, (a0+i0)
				; T(i+1,i) = 0;	Tcopy.pData[i*n + i-1]
				move r13.s, i0.s
				move r30.l, (a0+i0)
				; if fl == 0
				cmp #0, r6.s
					; q = i;
					j.ne QRSRAfterInnerLoop
					move.eq r7.s, r5.s
					j QRSREndInnerLoop
				j QRSRStartInnerLoop
			QRSRCond1False:
				; fl = 1;
				move #1, r6.s
				; p = i;
				move r7.s, r4.s

			QRSREndInnerLoop:
			dec r7.s, r7.s
			cmp #0, r7.s
			j.ne QRSRStartInnerLoop

; x 1 x 3 x x x 7 x 9 x 11 x 13 14 15 16 17 18 19 20 21 22 23 24 x 26 x x x x x
		; if q > 1
		QRSRAfterInnerLoop:
		cmp #1, r5.s
		j.le QRSREndMainLoop
			; [Qm,Rm] = qr_step(T(p:q,p:q));
			; size = q - p + 1;
			sub r4.s, r5.s, r8
			inc r8.s, r8.s ; size
			; формируем матрицу T(p:q,p:q)
			move g_Buffer8, a3.s ; T(p:q,p:q)
			; начальное смещение (p-1)* n + p - 1 = (p-1)*(n+1)
			dec r4.s, r6.s ; p - 1
			inc r0.s, r16.s ; n + 1
			mpuu r6.s, r16.s, r14.l ; (p-1)*(n+1)
			; инкремент адреса после внутреннего цикла n-q + p-1
			sub r5.s, r0.s, r9.s ; n-q
			add r9.s, r6.s, r9.s ; n-q + p-1
			; прибавляем к a0 начальное смещение
			move a0.s, r10.s
			add r10.s, r14.s, r10.s
			move r10.s, a0.s

			move #0, r12.s
			QRSRPartialCopyOfTOuterLoopBegin:
				do r8.s, QRSRPartialCopyOfTInnerLoopEnd
					; копируем элементы матрицы T
					move (a0)+, r10.l
				QRSRPartialCopyOfTInnerLoopEnd:
					move r10.l, (a3)+

				; увеличиваем a0 на n-q + p-1 (r9)
				move a0.s, r10.s
				add r10.s, r9.s, r10.s
				move r10.s, a0.s
				inc r12.s, r12.s
				cmp r12.s, r8.s
				j.ne QRSRPartialCopyOfTOuterLoopBegin

			; QRSimpleStep: ; здесь нет перехода, всё выполняется в теле функции
; r0.s - N
; r2.l - tol
; r4.s, r5.s - p, q
; r8.s - n
; x 1 x 3 x x 6 7 x 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 x 26 x x x x x
			dec r8.s, r1.s ; n - 1 ; здесь n - размер подматрицы T(p:q,p:q), записанный в r8
			dec r1.s, r7.s ; n - 2
			inc r8.s, r6.s ; n + 1
			mpuu r8.s, r8.s, r10.l ; n * n
			dec r10.s, r9.s ; n*n - 1
			dec r9.s, r9.s ; n*n - 2
; r1.s - n-1
; r7.s - n-2
; r6.s - n+1
; r10.s - n*n
; r9.s - n*n - 2
; x 1 x 3 x x x 7 x 9 x 11 12 13 14 15 16 17 18 19 20 21 22 23 24 x 26 x x x x x

			; d = (Tpq(n-1,n-1) - Tpq(n,n))/2;
			; // d = (Tcopy.pData[(n-2)*n + n-2] - Tcopy.pData[(n-1)*n + n-1]) / 2.0f;
			; // d = (Tcopy.pData[(n-2)*(n+1)] - Tcopy.pData[(n-1)*(n+1)]) / 2.0f;
			move g_Buffer8, a3.s
			mpuu r7.s, r6.s, r12.l ; (n-2)*(n+1)
			move r12.s, i3.s
			move (a3 + i3), r14.l ; Tpq(n-1,n-1)
			mpuu r1.s, r6.s, r16.l ; (n-1)*(n+1)
			move r16.s, i3.s
			move (a3 + i3), r18.l ; Tpq(n,n)
			fsub r18.l, r14.l, r14.l ; Tpq(n-1,n-1) - Tpq(n,n)
			fmpy r14.l, r31.l, r1.l ; d
; r18.l - Tpq(n,n)
; r1.l - d
; x x x 3 x x x 7 x 9 x 11 12 13 14 15 16 17 x 19 20 21 22 23 24 x 26 x x x x x

			; mu = Tpq(n,n) - Tpq(n,n-1)*Tpq(n,n-1)'/(d + d/norm(d)*sqrt(d*d' + Tpq(n,n-1)*Tpq(n,n-1)'));
			move r9.s, i3.s
			move (a3 + i3), r14.l ; Tpq(n, n-1) ; f1 = Tcopy.pData[n*n - 2];
			fmpy r14.l, r14.l, r3.l ; f1 = f1 * f1;
			; // mu = Tcopy.pData[n*n - 1] - f1/(d + d/f_abs(d)*sq_rt(d*d + f1));
			; d/f_abs(d) = +-1.0F
			andl 0x80000000, r1.l, r7.l
			move.eq 0x3F800000, r7.l ; +1.0
			move.ne 0xBF800000, r7.l ; -1.0
			inc r9.s, r9.s ; n*n + 1
			fmpy r1.l, r1.l, r9.l ; d*d
			fadd r9.l, r3.l, r9.l ; d*d + f1
			; sqrt(d*d + f1)
			; https://ru.wikipedia.org/wiki/Быстрый_обратный_квадратный_корень
			; y  = y * ( 1.5 - ( x/2.0 * y * y ) );
			fmpy r9.l, r31.l, r15.l ; x / 2.0
			finr r9.l, r9.l ; y
			do 5, QRSRSqrtLoop0End
				fmpy r15.l, r9.l, r11.l ; x/2 * y
				fmpy r11.l, r9.l, r11.l ; x/2 * y * y
				fsub r11.l, r29.l, r13.l ; 1.5 - ( x/2.0 * y * y )
			QRSRSqrtLoop0End:
				fmpy r9.l, r13.l, r9.l ; 1 / sqr(x) = y * ( 1.5 - ( x/2.0 * y * y ) );
			; теперь ищем 1/y
			; вычисляем второе приближение по формуле Ньютона-Рафсона y_{i+1} = y_i*(2 - x*y_i)
			fin r9.l, r11.l ; y
			do 5, QRSSEndSqrt
				fmpy r9.l, r11.l, r13.l ; x*y
				fsub r13.l, r27.l, r13.l ; 2 - x*y
			QRSSEndSqrt:
				fmpy r11.l, r13.l, r11.l ; y * (2 - x*y)
; r11.l - sqrt(d*d + f1)
			fmpy r7.l, r11.l, r11.l ; d/f_abs(d)*sq_rt(d*d + f1)
			fadd r1.l, r11.l, r11.l ; d + d/f_abs(d)*sq_rt(d*d + f1)
			; 1 / (d + d/f_abs(d)*sq_rt(d*d + f1))
			fin r11.l, r9.l
			do 5, QRSSEndInv
				fmpy r9.l, r11.l, r13.l
				fsub r13.l, r27.l, r13.l
			QRSSEndInv:
				fmpy r9.l, r13.l, r9.l
; r9.l - 1 / (d + d/f_abs(d)*sq_rt(d*d + f1))
			fmpy r9.l, r3.l, r9.l ; f1 / (d + d/f_abs(d)*sq_rt(d*d + f1))
			fsub r9.l, r18.l, r3.l ; mu
; r3.l - mu
; x x x x x x x 7 x 9 x 11 12 13 14 15 16 17 x 19 20 21 22 23 24 x 26 x x x x x

			; x = Tpq(1,1) - mu;
			move (a3), r12.l
			fsub r3.l, r12.l, r26.l
			; z = Tpq(2,1); // z = Tcopy.pData[n];
			move r8.s, i3.s
			move (a3 + i3), r12.l
; r26.l - x
; r12.l - z
; r18.l - free
; x x x x x x x 7 x 9 x 11 x 13 14 15 16 17 18 19 20 21 22 23 24 x x x x x x x

			; Qm = eye(n,n,'single');
			move g_Buffer2, a1.s
			move #0, r14.s
			QRSREyeOuterLoop1Begin:
				move r28.l, (a1)+
				;do r8.s, QRSREyeInnerLoop1End ; не работает при r8.s == 2
				;QRSREyeInnerLoop1End:
				;	move r30.l, (a1)+
				move #0, r15.s
				QRSREyeInnerLoop1Begin:
					move r30.l, (a1)+
					inc r15.s, r15.s
					cmp r15.s, r8.s
					j.ne QRSREyeInnerLoop1Begin
				inc r14.s, r14.s
				cmp r14.s, r1.s
				j.ne QRSREyeOuterLoop1Begin
			move r28.l, (a1)
			move g_Buffer2, a1.s

			; for k = 1:n-1 // for (k = 0; k < n-1; ++k)
; r18.s - k
			move #0, r18.s
			QRSRQRSimpleStepBeginMainLoop:
				; [c, s] = giv(x,z);

				; if abs(z) == 0
				andl 0x7fffffff, r12.l, r11.l ; |z|
				fsub r25.l, r11.l, r7.l ; |z| - FLOAT_EPS
				andl 0x80000000, r7.l, r7.l ; (|z| - FLOAT_EPS) < 0
				j.eq QRSRGivensAbsZNotZero
					move r28.l, r20.l ; c = 1
					move r30.l, r22.l ; s = 0
					j SQRSGivensEnd
; r20.l - c
; r22.l - s
; x x x x x x x 7 x 9 x 11 x 13 14 15 16 17 x 19 x 21 x 23 24 x x x x x x x
				; else (if abs(z) > 0)
				QRSRGivensAbsZNotZero:
					; if abs(z) > abs(x)
					andl 0x7fffffff, r26.l, r13.l ; |x|
					fsub  r13.l, r11.l, r11.l ; |z| - |x|
					andl 0x80000000, r11.l, r11.l ; (|z| - |x|) > 0
					j.ne SQRSGivensAbsAGreaterOrEqualThanAbsB
						; tau = -x/z
						fin r12.l, r11.l
						; вычисляем второе приближение по формуле Ньютона-Рафсона y_{i+1} = y_i*(2 - x*y_i)
						do 5, QRSRGivensInvZEnd
							fmpy r12.l, r11.l, r13.l ; x*y_i
							fsub r13.l, r27.l, r15.l ; 2 - x*y_i
						QRSRGivensInvZEnd:
							fmpy r11.l, r15.l, r11.l ; y_i*(2 - x*y_i)
						fmpy r26.l, r11.l, r11.l ; x/z
						fsub r11.l, r30.l, r11.l ; tau = -x/z
; r11.l - tau
						fmpy r11.l, r11.l, r13.l ; tau*tau
						fadd r13.l, r28.l, r13.l ; x = 1 + tau*tau'
						; s = 1/sqrt(1+tau*tau');
						fmpy r13.l, r31.l, r21.l ; / x/2
						finr r13.l, r22.l ; y
						do 5, QRSRSqrt1LoopEnd
							fmpy r21.l, r22.l, r15.l ; x/2 * y
							fmpy r15.l, r22.l, r15.l ; x/2 * y * y
							fsub r15.l, r29.l, r15.l ; 1.5 - ( x/2.0 * y * y )
						QRSRSqrt1LoopEnd:
							fmpy r15.l, r22.l, r22.l
						; c = s*tau;
						fmpy r22.l, r11.l, r20.l ; c
						j SQRSGivensEnd
					SQRSGivensAbsAGreaterOrEqualThanAbsB:
						; tau = -z/x
						fin r26.l, r11.l
						; вычисляем второе приближение по формуле Ньютона-Рафсона y_{i+1} = y_i*(2 - x*y_i)
						do 5, QRSRGivensInvXEnd
							fmpy r26.l, r11.l, r13.l ; x*y_i
							fsub r13.l, r27.l, r15.l ; 2 - x*y_i
						QRSRGivensInvXEnd:
							fmpy r11.l, r15.l, r11.l ; y_i*(2 - x*y_i)
						fmpy r12.l, r11.l, r11.l ; z/x
						fsub r11.l, r30.l, r11.l ; tau = -z/x
						; c = 1/sqrt(1+tau*tau');
						fmpy r11.l, r11.l, r13.l ; tau*tau
						fadd r13.l, r28.l, r13.l ; x = 1 + tau*tau'
						fmpy r13.l, r31.l, r21.l ; / x/2
						finr r13.l, r20.l ; y
						do 5, QRSRSqrt2LoopEnd
							fmpy r21.l, r20.l, r15.l ; x/2 * y
							fmpy r15.l, r20.l, r15.l ; x/2 * y * y
							fsub r15.l, r29.l, r15.l ; 1.5 - ( x/2.0 * y * y )
						QRSRSqrt2LoopEnd:
							fmpy r20.l, r15.l, r20.l ; c = 1 / sqr(x) = y * ( 1.5 - ( x/2.0 * y * y ) );
						; s = c*tau;
						fmpy r20.l, r11.l, r22.l ; s
				SQRSGivensEnd:
; x x x x x x x x x x x 11 x 13 14 15 16 17 x 19 20 21 22 23 24 x x x x x x x
				; G = eye(n,n,'single');
				move g_Buffer3d, a4.s
				move #0, r14.s
				QRSREyeOuterLoop2Begin:
					move r28.l, (a4)+
					;do r8.s, QRSREyeInnerLoop2End ; не работает при r8.s == 2
					;QRSREyeInnerLoop2End:
					;	move r30.l, (a4)+
					move #0, r15.s
					QRSREyeInnerLoop2Begin:
						move r30.l, (a4)+
						inc r15.s, r15.s
						cmp r15.s, r8.s
						j.ne QRSREyeInnerLoop2Begin
					inc r14.s, r14.s
					cmp r14.s, r1.s
					j.ne QRSREyeOuterLoop2Begin
				move r28.l, (a4)
				move g_Buffer3d, a4.s
				; G(k:k+1,k:k+1) = [c s; -s' c];
				; // G.pData[k*n + k] = c; // G.pData[k*(n+1)] = c;
				mpuu r18.s, r6.s, r14.l ; k*(n+1)
				move r14.s, i4.s
				move r20.l, (a4+i4)
				; // G.pData[k*n + k+1] = s; //  // G.pData[k*(n+1) + 1] = s;
				inc r14.s, r14.s
				move r14.s, i4.s
				move r22.l, (a4+i4)
				; // G.pData[(k+1)*n + k+1] = c; ; // G.pData[(k+1)*(n+1)] = c;
				inc r18.s, r19.s ; k+1
				mpuu r19.s, r6.s, r14.l
				move r14.s, i4.s
				move r20.l, (a4+i4)
				; // G.pData[(k+1)*n + k] = -s; // G.pData[(k+1)*n + k] = -s;
				dec r14.s, r14.s
				move r14.s, i4.s
				fsub r22.l, r30.l, r22.l ; -s
				move r22.l, (a4+i4)

; r20.l - free
; r22.l - free
; x x x x x x x x x x x 11 x 13 14 15 16 17 x 19 20 21 22 23 24 x x x x x x x
				; Tpq = G'*Tpq*G;
; G - g_Buffer3d, a4
; Tpq - g_Buffer8, a3
; G'*Tpq - g_Buffer4d, a5
				; G'*Tpq
				move g_Buffer3d, a4.s
				move g_Buffer8, a3.s
				move g_Buffer4d, a5.s
				move r8.s, i3.s ; инкремент адреса n
				move r8.s, i4.s ; инкремент адреса n
				move a3.s, r24.s ; начальное значение a3
				move a4.s, r25.s ; начальное значение a4
				move #0, r14.s ; счётчик внешнего цикла
				QRSRTGtTOuterLoopBegin:
					move #0, r15.s ; счётчик внутреннего цикла
					QRSRTGtTInnerLoopBegin:
						add r25.s, r14.s, r16.s
						move r16.s, a4.s
						add r24.s, r15.s, r16.s
						move r16.s, a3.s
						move #0, r16.l ; здесь накапливаем сумму
						do r8.s, QRSRTGtTAccLoopEnd
							move (a3)+i3, r20.l
							move (a4)+i4, r22.l
							fmpy r20.l, r22.l, r20.l
						QRSRTGtTAccLoopEnd:
							fadd r16.l, r20.l, r16.l
						move r16.l, (a5)+
						inc r15.s, r15.s
						cmp r15.s, r8.s
						j.ne QRSRTGtTInnerLoopBegin
					inc r14.s, r14.s
					cmp r14.s, r8.s
					j.ne QRSRTGtTOuterLoopBegin

				; Tpq = (G'*Tpq) * G
				move g_Buffer4d, a5.s
				move g_Buffer3d, a4.s
				move g_Buffer8, a3.s
				move a5.s, r24.s ; начальное значение a5
				move a4.s, r25.s ; начальное значение a4
				move #0, r14.s ; счётчик внешнего цикла
				QRSRTGtTGOuterLoopBegin:
					move #0, r15.s ; счётчик внутреннего цикла
					QRSRTGtTGInnerLoopBegin:
						mpuu r14.s, r8.s, r16.l
						add r24.s, r16.s, r16.s
						move r16.s, a5.s
						add r25.s, r15.s, r16.s
						move r16.s, a4.s
						move #0, r16.l ; здесь накапливаем сумму
						do r8.s, QRSRTGtTGAccLoopEnd
							move (a5)+, r20.l
							move (a4)+i4, r22.l
							fmpy r20.l, r22.l, r20.l
						QRSRTGtTGAccLoopEnd:
							fadd r16.l, r20.l, r16.l
						move r16.l, (a3)+
						inc r15.s, r15.s
						cmp r15.s, r8.s
						j.ne QRSRTGtTGInnerLoopBegin
					inc r14.s, r14.s
					cmp r14.s, r8.s
					j.ne QRSRTGtTGOuterLoopBegin

				; Qm = Qm*G;
				; сначала копируем Qm в g_Buffer4d
; x x x x x x x x x x x 11 x 13 14 15 16 17 x 19 20 21 22 23 24 x x x x x x x
; Qm - g_Buffer2, a1
; Qmcopy - g_Buffer4d, a5
; G - g_Buffer3d, a4
				move g_Buffer2, a1.s
				move g_Buffer4d, a5.s
				do r10.s, QRSRQCopyEnd
					move (a1)+, r12.l
				QRSRQCopyEnd:
					move r12.l, (a5)+

				move g_Buffer2, a1.s
				move g_Buffer3d, a4.s
				move g_Buffer4d, a5.s
				move a5.s, r24.s ; начальное значение a5
				move a4.s, r25.s ; начальное значение a4
				move #0, r14.s ; счётчик внешнего цикла
				QRSRTQGOuterLoopBegin:
					move #0, r15.s ; счётчик внутреннего цикла
					QRSRTQGInnerLoopBegin:
						mpuu r14.s, r8.s, r16.l
						add r24.s, r16.s, r16.s
						move r16.s, a5.s
						add r25.s, r15.s, r16.s
						move r16.s, a4.s
						move #0, r16.l ; здесь накапливаем сумму
						do r8.s, QRSRTQGAccLoopEnd
							move (a5)+, r20.l
							move (a4)+i4, r22.l
							fmpy r20.l, r22.l, r20.l
						QRSRTQGAccLoopEnd:
							fadd r16.l, r20.l, r16.l
						move r16.l, (a1)+
						inc r15.s, r15.s
						cmp r15.s, r8.s
						j.ne QRSRTQGInnerLoopBegin
					inc r14.s, r14.s
					cmp r14.s, r8.s
					j.ne QRSRTQGOuterLoopBegin

				; if k < n-1
				sub r1.s, r18.s, r24.s
				j.ge QRSRQRSimpleStepBeginMainLoop
					; x = T(k+1,k); // x = Tcopy.pData[(k+1)*n + k];
					move g_Buffer8, a3.s ; Tpq
					mpuu r19.s, r8.s, r14.l ; (k+1)*n
					add r14.s, r18.s, r14.s ; (k+1)*n + k
					move r14.s, i3.s
					move (a3+i3), r26.l
					; z = T(k+2,k);
					add r14.s, r8.s, r14.s ; (k+2)*n + k
					move r14.s, i3.s
					move (a3+i3), r12.l

				QRSRQRSimpleStepEndMainLoop:
				inc r18.s, r18.s
				cmp r18.s, r1.s
				j.ne QRSRQRSimpleStepBeginMainLoop:

; x x x x x x x x x x x 11 x 13 14 15 16 17 x 19 20 21 22 23 24 x x x x x x x
			; T(p:q,p:q) = Tpq;
			move g_Buffer1, a0.s ; T
			move g_Buffer8, a3.s ; Tpq
			dec r4.s, r14.s ; p - 1
			dec r5.s, r15.s ; q - 1
			sub r5.s, r0.s, r16.s
			add r16.s, r14.s, r16.s ; (n-q) + (p-1)
			sub r4.s, r5.s, r17.s
			inc r17.s, r17.s ; q - p + 1
			move a0.s, r20.s
			mpuu r14.s, r0.s, r22.l ; (p-1) * N
			add r20.s, r22.s, r20.s
			add r20.s, r14.s, r20.s
			move r20.s, a0.s ; T + (p-1)
			move #0, r20.s
			QRSRCopyRmToTLoopBegin:
				do r17.s, QRSRCopyRmToTInnerLoopEnd
					move (a3)+, r22.l
				QRSRCopyRmToTInnerLoopEnd:
					move r22.l, (a0)+
				move a0.s, r22.s
				add r22.s, r16.s, r22.s
				move r22.s, a0.s
				inc r20.s, r20.s
				cmp r20.s, r17.s
				j.ne QRSRCopyRmToTLoopBegin

			; Q = eye(N,N,'single');
			move g_Buffer3d, a5.s
			move #0, r20.s
			QRSREyeOuterLoop3Begin:
				move r28.l, (a5)+
				do r0.s, QRSREyeInnerLoop3End
				QRSREyeInnerLoop3End:
					move r30.l, (a5)+
				inc r20.s, r20.s
				cmp r20.s, r1.s
				j.ne QRSREyeOuterLoop3Begin
			move r28.l, (a5)

			; Q(p:q,p:q) = Qm;
			move g_Buffer2, a1.s ; Qm
			move g_Buffer3d, a5.s ; Q
			move a5.s, r20.s
			mpuu r14.s, r0.s, r22.l ; (p-1) * N
			add r20.s, r22.s, r20.s
			add r20.s, r14.s, r20.s
			move r20.s, a5.s ; Q + (p-1)
			move #0, r20.s
			QRSRCopyQmToQLoopBegin:
				do r17.s, QRSRCopyQmToQInnerLoopEnd
					move (a1)+, r22.l
				QRSRCopyQmToQInnerLoopEnd:
					move r22.l, (a5)+
				move a5.s, r22.s
				add r22.s, r16.s, r22.s
				move r22.s, a5.s
				inc r20.s, r20.s
				cmp r20.s, r17.s
				j.ne QRSRCopyQmToQLoopBegin

; x x x x x x x x x x x 11 x 13 14 15 16 17 x 19 20 21 22 23 24 x x x x x x x
			; U = U*Q;
			; копируем U в g_Buffer4d
			move g_Buffer7, a1.s ; U
			move g_Buffer4d, a4.s ; Ucopy
			mpuu r0.s, r0.s, r14.l ; N * N
			do r14.s, QRSRCopyULoopEnd
				move (a1)+, r16.l
			QRSRCopyULoopEnd:
				move r16.l, (a4)+

			move g_Buffer7, a1.s ; U
			move g_Buffer4d, a4.s ; UCopy
			move g_Buffer3d, a3.s ; Q
			move r0.s, i3.s
			move a4.s, r16.s ; начальное значение a4.s
			move a3.s, r17.s ; начальное значение a3.s
			move #0, r14.s
			QRSRQUQOuterLoopBegin:
				move #0, r15.s
				QRSRQUQInnerLoopBegin:
					move r16.s, a4.s
					add r17.s, r15.s, r19.s
					move r19.s, a3.s
					move #0, r20.l ; здесь накапливаем сумму
					do r0.s, QRSRQUQAccLoopEnd
						move (a4)+, r22.l
						move (a3)+i3, r24.l
						fmpy r22.l, r24.l, r11.l
					QRSRQUQAccLoopEnd:
						fadd r20.l, r11.l, r20.l
					move r20.l, (a1)+
					inc r15.s, r15.s
					cmp r15.s, r0.s
					j.ne QRSRQUQInnerLoopBegin
				add r16.s, r0.s, r16.s
				inc r14.s, r14.s
				cmp r14.s, r0.s
				j.ne QRSRQUQOuterLoopBegin

	QRSREndMainLoop:
	cmp #1, r5.s
	j.ne QRSRStartMainLoop

stop
