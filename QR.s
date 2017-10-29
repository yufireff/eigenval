 .global EigSymmTriag
 .global QRSimpleStep
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
 .global g_nColumns2
 .global g_ptr0

 .text
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Симметричный QR-алгоритм
; Источник: Голуб Дж., Ван Лоун Ч. Матричные вычисления: Пер. с англ. -
;              М.: Мир, 1999. - 548 с. (стр. 380 )
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;function [U,D] = QRsymm(T, tol)
;	n = size(T,1);
;	U = eye(n,n,'single');
;
;	p = 1;
;	q = n;
;	while q ~= 1
;		fl = 0;
;		for i = q-1:-1:1
;			if (abs(T(i+1,i)) <= tol*(abs(T(i,i)) + abs(T(i+1,i+1)))) % Cond1
;				T(i,i+1) = 0;
;				T(i+1,i) = 0;
;
;				if fl == 0
;					q = i;
;				else
;					break;
;				end
;			else
;				fl = 1;
;				p = i;
;			end
;		end
;
;		if q > 1
;			[Qm,Rm] = qr_step(T(p:q,p:q));
;			T(p:q,p:q) = Rm;
;
;			Q = eye(n, n);
;			Q(p:q,p:q) = Qm;
;			U = U*Q;
;		end
;	end
;	D = T;
;end
QRSymmReal:
; T находится в g_Buffer1
; n - в g_Size
; tol - в g_Factor
; на выходе получаем матрицы U (Q, g_Buffer2), D (S, g_Buffer5)
	move g_Size, a0.s
	move (a0), r0.l ; n
	move g_Factor, a0.s
	move (a0), r2.l ; tol
	move g_ptr0, a3.s
	move g_Buffer1, a0.s
	move g_Buffer2, a1.s
	move g_Buffer5, a2.s
	move g_p, a4.s ; p
	move g_q, a5.s ; q

	move #0, r30.l ; литерал 0.0f
	;	p = 1;
	move #1, r4.s ; p
	;	q = n;
	move r0.s, r5.s ; q
; занятые регистры ri.l
; x 1 x 3 x 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 x 31

	;	while q ~= 1
	QRSRStartMainLoop:
	;		fl = 0;
		move #1, r6.s ; f1
	;		for i = q-1:-1:1
		dec r5.s, r7.s ; i
; x 1 x 3 x 5 x 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 x 31
		QRSRStartInnerLoop:
	;			if (abs(T(i+1,i)) <= tol*(abs(T(i,i)) + abs(T(i+1,i+1))))
			inc r7.s, r8.s ; i+1
; x 1 x 3 x 5 x 7 x 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 x 31
			mpuu r7.s, r0.s, r10.l ; i*n
			mpuu r8.s, r0.s, r12.l ; (i+1)*n
; x 1 x 3 x 5 x 7 x 9 x 11 x 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 x 31
			add r12.s, r7.s, r11.s
			move r11.s, i0.s
			move (a0+i0), r14.l ; T(i+1,i)
			andl 0x7fffffff, r14.l, r1.l ; abs(T(i+1,i))
			add r10.s, r7.s, r11.s
			move r11.s, i0.s
			move (a0+i0), r14.l ; T(i,i)
			andl 0x7fffffff, r14.l, r3.l ; abs(T(i,i))
			add r12.s, r8.s, r11.s
			move r11.s, i0.s
			move (a0+i0), r14.l ; T(i+1,i+1)
			andl 0x7fffffff, r14.l, r5.l ; abs(T(i+1,i+1))
			fadd r3.l, r5.l, r7.l ; abs(T(i,i)) + abs(T(i+1,i+1))
			fmpy r2.l, r7.l, r9.l ; tol*(abs(T(i,i)) + abs(T(i+1,i+1)))
			fsub r1.l, r9.l, r11.l ; abs(T(i+1,i)) - tol*(abs(T(i,i)) + abs(T(i+1,i+1)))
			; проверяем знак получившегося числа. если старший бит = 1, число отрицательное
			cmpl 0x7fffffff, r11.l ; если результат > 0, то разность отрицательна, и изначальное условие false
			j.gt QRSRCond1False
	;				T(i,i+1) = 0;
				add r10.s, r8.s, r11.s
				move r11.s, i0.s
				move (a0+i0), r30.l
	;				T(i+1,i) = 0;
				add r12.s, r7.s, r11.s
				move r11.s, i0.s
				move (a0+i0), r30.l
	;				if fl == 0
				cmp #0, r6.s
	;					q = i;
					move.eq r7.s, r5.s
					j.ne QRSRAfterInnerLoop
			QRSRCond1False:
	;				fl = 1;
				move #1, r6.s
	;				p = i;
				move r7.s, r4.s
			dec r7.s, r7.s
			cmp #1, r7.s
			j.ne QRSRStartInnerLoop

	;		if q > 1
	QRSRAfterInnerLoop:
		cmp #1, r5.s
		j.gt QRSREndMainLoop
		move r4.s, (a4)
		move r5.s, (a5)
	;			[Qm,Rm] = qr_step(T(p:q,p:q));
		bs QRSimpleStep
; x 1 x 3 x 5 x 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 x 31

	;			T(p:q,p:q) = Rm;
		move g_Buffer7, a4.s
		move a4.s, r13.s
		move r4.s, r8.s ; p
		sub r5.s, r4.s, r12.s ; q - p
		QRSRStartCopyRmLoop:
			mpuu r8.s, r0.s, r10.l ; p * n
			add r4.s, r10.s ; p*n  + p
			add r13.s, r10.s
			move r10.s a0.s
			do r12.s, QRSREndCopyRmInnerLoop
				move (a4)+, r14.l
			QRSREndCopyRmInnerLoop:
				move r14.l, (a0)+
			inc r8.s, r8.s
			cmp r8.s, r5.s
			j.ne QRSRStartCopyRmLoop

	;			Q = eye(n, n);
		move (a3), a1.s
		bs MatrixRealEye
	;			Q(p:q,p:q) = Qm;
	;			U = U*Q;

	QRSREndMainLoop:
		cmp #1, r5.s
		j.ne QRSRStartMainLoop

	rts

; Неявный симметричный QR - шаг неприводимой
; трехдиагональной матрицы со сдвигом Уилкинсона
;function [Q, T] = qr_step(T)
;	n = size(T,2);
;	d = (T(n-1,n-1) - T(n,n))/2;
;
;	mu = T(n,n) - T(n,n-1)*T(n,n-1)'/(d + d/norm(d)*sqrt(d*d' + T(n,n-1)*T(n,n-1)'));
;
;	x = T(1,1) - mu;
;	z = T(2,1);
;	Q = eye(n,n,'single');
;	for k = 1:n-1
;		[c, s] = giv(x,z);
;		G = eye(n,n,'single');
;		G(k:k+1,k:k+1) = [c s; -s' c];
;		T = G'*T*G;
;		Q = Q*G;
;
;		if k < n-1
;			x = T(k+1,k);
;			z = T(k+2,k);
;		end
;	end
;end
QRSimpleStep:
; T находится в А0
	move g_Size, a7.s
	move (a7), r0.l ; N
	dec r0.s, r1.s ; N-1

	;	d = (T(n-1,n-1) - T(n,n))/2;
	mpuu r1.s, r1.s, r2.l
	;move r2.s,
