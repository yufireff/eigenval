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

	; ���������� ��������� ������
	move a2, r5.l
	move a3, r6.l

	move #0, r8.l; eye
	move 0x3f800000, r10.l; 1.0f

	; ������ �������� r0, r2, r4, r5, r6, r8, r10

	move #0, r1.l; ������� �����
	; ���������� ������������� �������
	CMMStartAssociatedMatrixRowsLoop:
		move #0, r3.l; ������� ��������
		CMMStartAssociatedMatrixColumnsLoop1: ; ����� ������� � �������
			move (a0)+, r14.l
			move r14.l, (a2)+
			move (a1)+, r14.l
			move r14.l, (a3)+
			incl r3.l, r3.l
			cmpl r3.l, r0.l
			j.ne CMMStartAssociatedMatrixColumnsLoop1
			move #0, r3.l
		CMMStartAssociatedMatrixColumnsLoop2: ; ����� ��������� �������
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

	; ���������� ������ �� ��������� ��������
	move r5.l, a2
	move r6.l, a3

	; �������� �������������� ������� � ������������������ ����
	move #0, r1.l; ������� �����
	subl #2, r0.l, r25.l ; r25 = N - 2
	CMMStartLowerLeftMainLoop:
		; ���������� ������������ �������
		mpuu r1.s, r4.s, r3.l ; r3 = 2N * i
		nop ; ��������� ���������� � 2 �����
		addl r1.l, r3.l, r7.l ; r7 = 2N*i + i
		move r3.s, i2.s
		move r3.s, i3.s
		move (a2 + i2), r12.l ; re diag
		move (a3 + i3), r14.l ; im diag
		; ������ ������������ ������� ������ 1
		; ��������� ����� ������������� ��������
		fmpy r12.l, r12.l, r13.l ; re^2 diag
		fmpy r14.l, r14.l, r15.l ; im^2 diag. ���������� ����� �����, ����� �� ��������� nop
		move #0, r9.l ; �������� ������ �����
		; ���������� ������ ������ ������
		addl r5.l, r3.l, r11.l
		;nop
		move r11.l, a2
		addl r6.l, r3.l, r27.l
		fadd r15.l, r13.l ; r13 = re^2 diag + im^2 diag. �������� ���������� ���� � ����� �����������
		move r27.l, a3.l
		move #1, i2.s
		move #1, i3.s
		fin r13.l, r13.l ; r13 = 1 / (re^2 diag + im^2 diag). ������� ���������� ����� � ����� �����������
		CMMStartDiagonalOneLoop:
			; ����� ������� ������� �� ������������
			move (a2), r16.l ; re
			move (a3), r18.l ; im
			fmpy r12.l, r16.l, r15.l ; re_i * re_diag
			fmpy r14.l, r18.l, r17.l ; im_i * im_diag
			;nop
			fadd r17.l, r15.l ; r15 = re_i * re_diag + im_i * im_diag
			; ����� ���� ������ nop, ������� ������ ������������
			fmpy r12.l, r18.l, r19.l ; re diag * im i
			fmpy r16.l, r14.l, r21.l ; re i * im diag
			fmpy r15.l, r13.l, r20.l ; ����� �������������� ����� �� �����
			;fsub r19.l, r21.l, r23.l ; ��������� ������ �����
			fsub r21.l, r19.l, r23.l ; ��������� ������ �����
			move r20.l, (a2)+ ; ���� ���������� ��������� ������ �����, ���������� ��������������
			fmpy r23.l, r13.l, r20.l ; ����� ������ ����� �� �����
			incl r9.l, r9.l ; ���� ������� ������ �����, ����������� ������� �����
			move r20.l, (a3)+ ; ���������� ������ ����� � �������
			cmpl r9.l, r4.l ; ��������� ������� ������ �� �����
			j.ne CMMStartDiagonalOneLoop

		; ���������, �� ��������� �� ��� ������. ���� ��, ������� ��������
		cmpl r1.l, r25.l
		j.eq CMMEndLowerLeftMainLoop

		; �������� ������� �������� �������
		; ��� ����� �������� �� ���� ��������� ����� ������� � �������������
		; ����� 2 ������:
		; �0, �1 - ��� ���������� ������. ���������� � base + 2N * i
		; a2, a3 - ��� ������� ����������. ���������� � base + 2N * j, j = [i+1, N]
		; ������ ���������� ������ ��������� � ��������� r11, r27
		; �������� � ������� ���������� ��������
		addl r1.l, r11.l
		addl r1.l, r27.l
		move #0, r30.l
		inc r1.s, r30.s ; �������� ���� � i + 1; ���� r1
		CMMStartZeroingColumnOuterLoop1:
			move r11.l, a0.l
			move r27.l, a1.l
			; ����� �2 � �3
			mpuu r4, r30, r14
			nop
			addl r1.l, r14.l
			nop
			addl r5.l, r14.l, r15.l
			addl r6.l, r14.l, r16.l
			move r15.l, a2.l
			move r16.l, a3.l
			; ���������� ������� �������, �� ������� ���������� ������
			move (a2), r14.l ; re 1
			move (a3), r18.l ; im 1
			; ���� �� ������
			move r30.l, r17.l
			CMMStartZeroingColumnInnerLoop1:
				; ������� ���������� ������
				move (a0)+, r16.l ; re 2
				move (a1)+, r26.l ; im 2
				; ������� ������� ������
				move (a2), r28.l ; re 3
				move (a3), r24.l ; im 3
				; ������������ ����������� �������� �� �������
				fmpy r14.l, r16.l, r19.l ; re1 * re2
				fmpy r18.l, r26.l, r20.l ; im1 * im2
				fmpy r14.l, r26.l, r21.l ; re1 * im2
				fmpy r18.l, r16.l, r22.l ; re2 * im1
				fsub r20.l, r19.l, r20.l; re1 * re2 - im1 * im2
				fadd r21.l, r22.l, r22.l; re1 * im2 + re2 * im1
				; �������� ��� �� �������� ��������
				fsub r20.l, r28.l, r20.l
				fsub r22.l, r24.l, r22.l
				; ���������� ���������
				move r20.l, (a2)+
				move r22.l, (a3)+

				incl r17.l, r17.l ; ����� �� ����� ��������� ������
				cmpl r17.l, r4.l
				j.ne CMMStartZeroingColumnInnerLoop1

			incl r30.l, r30.l ; ����� �� ����� ��������� �������
			cmpl r30.l, r0.l
			j.ne CMMStartZeroingColumnOuterLoop1

		incl r1.l, r1.l ; ����� �� ����� ��������� ������� ������������
		cmpl r1.l, r0.l
		j.ne CMMStartLowerLeftMainLoop
	CMMEndLowerLeftMainLoop:

	stop
