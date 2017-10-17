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

; ��������� �������������� ������ A' * B * A
RealMatrixMultAtBA:
	move g_pReal1, a0.s
	move g_pReal2, a1.s
	move g_pReal3, a2.s
	move g_pImag1, a3.s ; ���������� ����� Imag1 ��� ������������� ����������
	; ��������� ������� �������
	move a0.s, r0.s
	move a1.s, r1.s
	move a2.s, r2.s
	move a3.s, r3.s

	move g_Size, a4.l
    move (a4), r4.l
; ������� �������� ri.l
; x 1 x 3 x 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

	; ������� ������� ������������ A' * B
	; ����� ��� ������
	; ��� �0 � �1 - ������ �������
	; ��� �3 - 1
	move r4.s, i0.s
	move r4.s, i1.s
	move #1, i3.s

	move #0, r7.s ; i, ������� ����� � �������������� �������
	RMMAtBAStartRowsLoop1: ; � ����� �������� ��� ����� ��������
		move #0, r9.s ; j, ������� �������� � �������������� �������
		RMMAtBAStartColumnsLoop1:
			; ��������� ��������� ������ a0 = base + i
			add r0.s, r7.s, r12.s
			move r12.s, a0.s
			; ��������� ��������� ������ a1 = base + j
			add r1.s, r9.s, r12.s
			move r12.s, a1.s
			move #0, r18.l ; ��������� �����
			do r4.s, RMMAtBAEndInnerLoop1
				move (a0)+i0, r12.l
				move (a1)+i1, r14.l
				fmpy r12.l, r14.l, r16.l
			RMMAtBAEndInnerLoop1:
				fadd r16.l, r18.l
				move r18.l, (a3)+
			inc r9.s, r9.s ; ������� ������ �� ����� �� ��������
			cmp r9.s, r4.s
			j.ne RMMAtBAStartColumnsLoop1
		inc r7.s, r7.s ; ������� ������ �� ����� �� �������
		cmp r7.s, r4.s
		j.ne RMMAtBAStartRowsLoop1

	; ������ ��������� �������� �� �
	move #1, i2.s
	move #0, r5.s ; i, ������� ����� � �������������� �������
	RMMAtBAStartRowsLoop2:
		move #0, r6.s ; j, ������� �������� � �������������� �������
		RMMAtBAStartColumnsLoop2:
			; ��������� ��������� ������ a3 = base + i* N
			mpuu r5.s, r4.s, r8.l
			add r3.s, r8.s
			move r8.s, a3.s
			; ��������� ��������� ������ a0 = base + j
			add r0.s, r6.s, r8.s
			move r8.s, a0.s
			; ��������� �����
			move #0, r10.l
; x 1 x 3 x 5 x 7 x 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
			do r4, RMMAtBAEndInnerLoop2
				move (a3)+, r12.l
				move (a0)+i0, r14.l
				fmpy r12.l, r14.l, r3.l
			RMMAtBAEndInnerLoop2:
				fadd r3.l, r10.l
			; ��������� ���������
			move r10.l, (a2)+

			inc r6.s, r6.s ; ������� ������ �� ����� �� ��������
			cmp r6.s, r4.s
			j.ne RMMAtBAStartColumnsLoop2
		inc r5.s, r5.s ; ������� ������ �� ����� �� �������
		cmp r5.s, r4.s
		j.ne RMMAtBAStartRowsLoop2
	stop

