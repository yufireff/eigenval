 .global MatrixRealEye

 .text
	MatrixRealEye:
	move g_ptr0, a1.s
	move (a1.s), r0.s
	move r0.s, a0.s
	move g_Size, a7.s
	move (a7), r0.l ; N
	dec r0.s, r1.s ; N-1
	move #0, r2.s
	move 0x3f800000, r4.l
	move #0, r6.l
	MREOuterLoopBegin:
		move r4.l, (a0)+
		do r0.s, MREInnerLoopEnd
		MREInnerLoopEnd:
			move r6.l, (a0)+
		inc r2.s, r2.s
		cmp r2.s, r1.s
		j.ne MREOuterLoopBegin
	move r4.l, (a0)
	rts

