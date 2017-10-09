.global ComplexMatrixInv
.global g_pRealDoubleSize
.global g_pImagDoubleSize

.text

ComplexMatrixInv:
	move #0, a0.l
	move g_nRows1, a0
	move (a0), r0.l ; r0.s = N
	mpuu #2, r0.s, r2.l
	move r2.s, r1.s ; r1.s = 2N
	; ������ ���������� �������� ����� � r0.l, ��� ��� r0.l = N | (2N << 16)
; ������� �������� ri.l
; x 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	move #0, a0.l
	move g_pReal1, a0
	move #0, a1.l
	move g_pImag1, a1
	move #0, a2.l
	move g_pReal2, a2
	move #0, a3.l
	move g_pImag2, a3

	; ���������� ��������� ������
	move a2, r2.s ; r2.s = a2
	move a3, r3.s ; r3.s = a3
; x 1 x 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	move #0, r4.l ; ������� 0.0f
	move 0x3f800000, r6.l ; ������� 1.0f
; x 1 x 3 x 5 x 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	; ���������� ������������� �������
	move #0, r8.s ; ������� �����
	CMMStartAssociatedMatrixRowsLoop:
		move #0, r9.s ; ������� ��������
; x 1 x 3 x 5 x 7 x 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
		CMMStartAssociatedMatrixColumnsLoop1: ; ����� ������� � �������
			move (a0)+, r10.l
			move r10.l (a2)+
			move (a1)+, r10.l
			move r10.l (a1)+
			inc r9.s, r9.s
			cmp r9.s, r0.s
			j.ne CMMStartAssociatedMatrixColumnsLoop1

		move #0, r9.s
		CMMStartAssociatedMatrixColumnsLoop2: ; ����� ��������� �������
			cmp r8.s, r9.s
			move.eq r4.l, (a2)+
			move.ne r6.l, (a2)+
			move r4.l, (a3)+
			inc r9.s, r9.s
			cmp r9.s, r0.s
			j.ne CMMStartAssociatedMatrixColumnsLoop2

		inc r8.s, r8.s
		cmp r8.s, r0.s
		j.ne CMMStartAssociatedMatrixRowsLoop
; x 1 x 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	; ���������� ������ �� ��������� ��������
	move r2.s, a2
	move r3.s, a3

	; �������� �������������� ������� � ������������������ ����
	move #0, r8.s ; ������� ����� i
	sub #2, r0.s, r4.s ; r4.s = N - 2
	move 0x40000000, r5.l ; ������� 2.0f
	CMMStartLowerLeftMainLoop:
		; ���������� ������������ �������
		mpuu r8.s, r1.s, r4.l ; r4.s = 2N * i
		add r4.s, r8.s, r6.s ; r6.s = 2N+i + i
; x 1 x 3 x x x 7 x 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
		move r4.s, i2
		move r4.s, i3
		move (a2 + i2), r10.l ; r10 = re_diag
		move (a3 + i3), r12.l ; r12 = im_diag

		; ��������� ����� ������������� ��������
		fmpy r10.l, r10.l, r1.l
		fmpy r12.l, r12.l, r3.l
		fadd r3.l, r1.l ; r1 = re^2 diag + im^2 diag
; x x x 3 x x x 7 x 9 x 11 x 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
		fin r1.l r3.l ; r1.l - ������ ����������� 1 / norm(diag)
		; ��������� ������ ����������� �� ������� �������-������� y_{i+1} = y_i*(2 - x*y_i)
		; ���� �������� �� ������, ����� ����� ��������� �������� � �����
		fmpy r1.l, r3.l, r7.l
		fsub r7.l, r5.l, r7.l
		fmpy r7.l, r3.l
; x x x x x x x 7 x 9 x 11 x 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
		; ���������� ������ ������ ������ (start + 2N*i)
		add r2.s, r4.s, r14.s
		add r3.s, r4.s, r15.s
; x x x x x x x 7 x 9 x 11 x 13 x 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
		move #1, i2
		move #1, i3

		; ������ ������������ ������� ������ 1
		move #0, r16.s
; x x x x x x x 7 x 9 x 11 x 13 x 15 x 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
		CMMStartDiagonalOneLoop:
			; ����� ������� ������� �� ������������
			move (a2), r18.l ; re_i
			move (a3), r20.l ; im_i
; x x x x x x x 7 x 9 x 11 x 13 x 15 x 17 x 19 x 21 22 23 24 25 26 27 28 29 30 31
			fmpy r18.l, r10.l, r7.l ; re i * re diag
			fmpy r20.l, r12.l, r9.l ; im_i * im_diag
			fadd r9.l, r7.l ; r7 = re_i * re_diag - im_i * im_diag
			fmpy r7.l, r3.l, r22.l; ����� �������������� ����� �� �����
			move r22.l, (a2)+ ; ���������� ���������

			fmpy r18.l, r12.l, r7.l ; re i * im diag
			fmpy r20.l, r10.l, r9.l ; im_i * re_diag
			fadd r9.l, r7.l ; r7 = re_i * im_diag - im_i * re_diag
			fmpy r7.l, r3.l, r22.l; ����� ������ ����� �� �����
			move r22.l, (a3)+ ; ���������� ���������

			inc r16.s, r16.s
			cmp r16.s, r1.s
			j.ne CMMStartDiagonalOneLoop

		; ���������, �� ��������� �� ��� ������. ���� ��, ������� ��������
		cmp r1.s, r4.s
		j.eq CMMEndLowerLeftMainLoop

		; �������� ������� �������� �������
		; ��� ����� �������� �� ���� ��������� ����� ������� � �������������
		; ����� 2 ������:
		; �0, �1 - ��� ���������� ������. ���������� � base + 2N * i
		; a2, a3 - ��� ������� ����������. ���������� � base + 2N * j, j = [i+1, N]
		; ������ ���������� ������ ��������� � ��������� r14, r15
		; �������� � ������� ���������� ��������
		add r8.s, r14.s
		add r8.s, r15.s
		inc r8.s, r22.s ; �������� ���� � i + 1
; x x x x x x x 7 x 9 x 11 x 13 x 15 x 17 x 19 x 21 x 23 24 25 26 27 28 29 30 31
		CMMStartZeroingColumnOuterLoop1:
			move r14.s, a0
			move r15.s, a1
			; ����� �2 � �3
			add r4.s, r0.s, r24.s ; 2N*(i+1)
			add r24.s, r2.s, r25.s
			move r25.s a2
			add r24.s, r3.s, r25.s
			move r25.s, a3
			; ���������� ������� �������, �� ������� ���������� ������
			move (a2), r24.l ; ������� move ��� ������ � �������� �������� - ������ ������...
			move r24.l, r7.l ; re 1
			move (a3), r24.l 
			move r24.l, r9.l ; im 1
; x x x x x x x x x x x 11 x 13 x 15 x 17 x 19 x 21 x 23 24 25 x 27 28 29 30 31
			; ���� �� ������
			move r22.s, r24.s
; x x x x x x x x x x x 11 x 13 x 15 x 17 x 19 x 21 x 23 x 25 x 27 28 29 30 31
			CMMStartZeroingColumnInnerLoop1:
				; ������� ���������� ������
				move (a0)+, r26.l ; �� ������� ������ ���������, ������� ���� �� ����� ������ � 2 �����
				move r26.l, r11.l ; re 2
				move (a1)+, r30.l ; im 2
				; ������� ������� ������
				move (a2), r26.l ; re 3
				move (a3), r28.l ; im 3 
; x x x x x x x x x x x x x 13 x 15 x 17 x 19 x 21 x 23 x 25 x 27 x 29 x 31
				; ������������ ����������� �������� �� �������
				fmpy r7.l, r11.l, r19.l ; re1 * re2
				fmpy r9.l, r26.l, r21.l ; im1 * im2
				fmpy r7.l, r30.l, r23.l ; re1 * im2
				fmpy r11.l, r9.l, r25.l ; re2 * im1
				fsub r21.l, r19.l ; re1 * re2 - im1 * im2
				fadd r25.l, r23.l ; re1 * im2 + re2 * im1
				; �������� ��� �� �������� ��������
				fsub r19.l, r26.l
				fsub r23.l, r28.l
				; ���������� ���������
				move r26.l, (a2)+
				move r28.l, (a3)+

				inc r24.s, r24.s ; ����� �� ����� ��������� ������
				cmp r24.s, r1.s
				j.ne CMMStartZeroingColumnInnerLoop1
			inc r22.s, r22.s ; ����� �� ����� ��������� �������
			cmp r22.s, r0.s
			j.ne CMMStartZeroingColumnOuterLoop1
		inc r8.s, r8.s
		cmp r8.s, r0.s
		j.ne CMMStartLowerLeftMainLoop
	CMMEndLowerLeftMainLoop:

	stop